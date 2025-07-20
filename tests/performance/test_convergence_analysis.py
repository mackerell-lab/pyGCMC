#!/usr/bin/env python3
"""
收敛性分析 - 对比FBP和SCF的收敛行为
"""

import pygcmc
import numpy as np
import time

def create_stable_water_system(n_waters):
    """创建相对稳定的水系统"""
    # 使用较大间距确保初始配置合理
    spacing = 0.5  # nm
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

def analyze_convergence(state, n_waters, algorithm, max_iterations=1000):
    """分析算法的收敛行为"""
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
    
    # 测试不同容差
    tolerances = [100.0, 10.0, 1.0, 0.1]
    results = []
    
    for tol in tolerances:
        # 设置参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = max_iterations
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algorithm == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置Drude位置
        state_copy = state.copy()
        for i in range(n_waters):
            state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
            state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
            state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
        
        # 计时和运行
        start = time.time()
        try:
            energy = force.calculateEnergySCF(state_copy)
            elapsed = (time.time() - start) * 1000
            converged = True
        except Exception as e:
            energy = float('nan')
            elapsed = (time.time() - start) * 1000
            converged = False
        
        # 计算最终位移
        max_disp = 0
        avg_disp = 0
        for i in range(n_waters):
            dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
            dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
            dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
            disp = np.sqrt(dx*dx + dy*dy + dz*dz)
            max_disp = max(max_disp, disp)
            avg_disp += disp
        avg_disp /= n_waters
        
        results.append({
            'tolerance': tol,
            'time': elapsed,
            'energy': energy,
            'converged': converged,
            'avg_disp': avg_disp * 1000,  # pm
            'max_disp': max_disp * 1000   # pm
        })
    
    return results

def compare_iteration_efficiency():
    """比较每次迭代的效率"""
    print("迭代效率对比测试")
    print("="*80)
    
    n_waters = 64  # 中等大小系统
    state = create_stable_water_system(n_waters)
    
    # 固定迭代次数，测试单次迭代时间
    fixed_iterations = [10, 50, 100]
    
    print(f"\n{n_waters}水分子系统，测试固定迭代次数的时间")
    print("-"*60)
    print(f"{'迭代次数':<10} {'SCF时间(ms)':<15} {'FBP时间(ms)':<15} {'每迭代SCF(ms)':<15} {'每迭代FBP(ms)':<15}")
    print("-"*60)
    
    for n_iter in fixed_iterations:
        # 创建force对象
        force_scf = pygcmc.DrudeForce()
        force_fbp = pygcmc.DrudeForce()
        
        # 添加Drude粒子
        charge = -1.71636
        k_spring = 418400.0
        polarizability = 1.71636**2 * 138.935456 / k_spring
        
        for force in [force_scf, force_fbp]:
            for i in range(n_waters):
                force.addParticle(
                    drudeIndex=5*i+1, parentIndex=5*i,
                    aniso1Index=-1, aniso2Index=-1,
                    aniso3Index=-1, aniso4Index=-1,
                    charge=charge, polarizability=polarizability,
                    aniso12=1.0, aniso34=1.0
                )
            
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    force.addScreenedPair(i, j, 1.3)
        
        # 设置参数 - 使用很大的容差确保不会提前收敛
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e10  # 永远不会收敛
        params.maxIterations = n_iter
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        
        force_scf.setSCFParameters(params)
        force_scf.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        force_fbp.setSCFParameters(params)
        force_fbp.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 测试SCF
        state_scf = state.copy()
        for i in range(n_waters):
            state_scf.atoms[5*i+1].x = state_scf.atoms[5*i].x
            state_scf.atoms[5*i+1].y = state_scf.atoms[5*i].y
            state_scf.atoms[5*i+1].z = state_scf.atoms[5*i].z
        
        start = time.time()
        force_scf.calculateEnergySCF(state_scf)
        scf_time = (time.time() - start) * 1000
        
        # 测试FBP
        state_fbp = state.copy()
        for i in range(n_waters):
            state_fbp.atoms[5*i+1].x = state_fbp.atoms[5*i].x
            state_fbp.atoms[5*i+1].y = state_fbp.atoms[5*i].y
            state_fbp.atoms[5*i+1].z = state_fbp.atoms[5*i].z
        
        start = time.time()
        force_fbp.calculateEnergySCF(state_fbp)
        fbp_time = (time.time() - start) * 1000
        
        print(f"{n_iter:<10} {scf_time:<15.1f} {fbp_time:<15.1f} {scf_time/n_iter:<15.2f} {fbp_time/n_iter:<15.2f}")

def main():
    """主测试函数"""
    print("FBP vs SCF 收敛性分析")
    print("="*80)
    
    # 测试128水系统
    n_waters = 128
    state = create_stable_water_system(n_waters)
    
    print(f"\n分析{n_waters}水分子系统的收敛行为")
    print("-"*80)
    
    # 分析SCF
    print("\nSCF收敛性:")
    print(f"{'容差':<10} {'时间(ms)':<12} {'能量(kJ/mol)':<15} {'平均位移(pm)':<12} {'最大位移(pm)':<12}")
    print("-"*70)
    
    scf_results = analyze_convergence(state, n_waters, "SCF", max_iterations=500)
    for r in scf_results:
        if r['converged']:
            print(f"{r['tolerance']:<10.1f} {r['time']:<12.1f} {r['energy']:<15.2f} {r['avg_disp']:<12.2f} {r['max_disp']:<12.2f}")
        else:
            print(f"{r['tolerance']:<10.1f} {r['time']:<12.1f} {'未收敛':<15} {r['avg_disp']:<12.2f} {r['max_disp']:<12.2f}")
    
    # 分析FBP
    print("\nFBP收敛性:")
    print(f"{'容差':<10} {'时间(ms)':<12} {'能量(kJ/mol)':<15} {'平均位移(pm)':<12} {'最大位移(pm)':<12}")
    print("-"*70)
    
    fbp_results = analyze_convergence(state, n_waters, "FBP", max_iterations=500)
    for r in fbp_results:
        if r['converged']:
            print(f"{r['tolerance']:<10.1f} {r['time']:<12.1f} {r['energy']:<15.2f} {r['avg_disp']:<12.2f} {r['max_disp']:<12.2f}")
        else:
            print(f"{r['tolerance']:<10.1f} {r['time']:<12.1f} {'未收敛':<15} {r['avg_disp']:<12.2f} {r['max_disp']:<12.2f}")
    
    # 对比分析
    print("\n\n收敛性对比分析:")
    print("-"*60)
    
    for i, tol in enumerate([100.0, 10.0, 1.0, 0.1]):
        scf = scf_results[i]
        fbp = fbp_results[i]
        
        if scf['converged'] and fbp['converged']:
            speedup = scf['time'] / fbp['time']
            energy_diff = abs(fbp['energy'] - scf['energy'])
            
            print(f"\n容差 {tol}:")
            print(f"  时间: SCF {scf['time']:.1f} ms, FBP {fbp['time']:.1f} ms")
            print(f"  加速: {speedup:.2f}x")
            print(f"  能量差: {energy_diff:.2f} kJ/mol")
            print(f"  位移: SCF {scf['avg_disp']:.2f} pm, FBP {fbp['avg_disp']:.2f} pm")
    
    # 测试迭代效率
    print("\n")
    compare_iteration_efficiency()

if __name__ == "__main__":
    main()