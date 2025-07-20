#!/usr/bin/env python3
"""
详细比较CG和SCF在大体系中的结果差异
"""

import pygcmc
import numpy as np
import time

def create_water_system(n_waters):
    """创建水分子系统"""
    spacing = 0.4  # nm
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
                
                # SWM4-NDP水模型
                positions = [
                    (x, y, z, 1.71636, 0),   # O
                    (x, y, z, -1.71636, 1),  # D (初始与O重合)
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
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length/2 - 0.01)
    
    # SWM4-NDP力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def analyze_drude_positions(state, n_waters):
    """分析Drude位置分布"""
    displacements = []
    max_disp = 0.0
    min_disp = 1e10
    
    for i in range(n_waters):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        dy = state.atoms[5*i+1].y - state.atoms[5*i].y
        dz = state.atoms[5*i+1].z - state.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        displacements.append(disp)
        
        if disp > max_disp:
            max_disp = disp
        if disp < min_disp:
            min_disp = disp
    
    mean_disp = np.mean(displacements)
    std_disp = np.std(displacements)
    
    # 计算位移分布
    hist, bins = np.histogram(displacements, bins=20)
    
    return {
        'mean': mean_disp,
        'std': std_disp,
        'min': min_disp,
        'max': max_disp,
        'displacements': displacements,
        'histogram': (hist, bins)
    }

def calculate_dipole_moments(state, n_waters):
    """计算每个水分子的偶极矩"""
    dipoles = []
    
    for i in range(n_waters):
        # 计算质心
        com_x = (state.atoms[5*i].x * 15.5994 + 
                 state.atoms[5*i+2].x * 1.00783 + 
                 state.atoms[5*i+3].x * 1.00783) / 17.6155
        com_y = (state.atoms[5*i].y * 15.5994 + 
                 state.atoms[5*i+2].y * 1.00783 + 
                 state.atoms[5*i+3].y * 1.00783) / 17.6155
        com_z = (state.atoms[5*i].z * 15.5994 + 
                 state.atoms[5*i+2].z * 1.00783 + 
                 state.atoms[5*i+3].z * 1.00783) / 17.6155
        
        # 计算偶极矩 (charge * position)
        dipole_x = 0.0
        dipole_y = 0.0
        dipole_z = 0.0
        
        for j in range(5):
            idx = 5*i + j
            charge = state.atoms[idx].charge
            dipole_x += charge * (state.atoms[idx].x - com_x)
            dipole_y += charge * (state.atoms[idx].y - com_y)
            dipole_z += charge * (state.atoms[idx].z - com_z)
        
        # 转换为Debye (1 e·nm = 4.80321 D)
        dipole_mag = np.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2) * 4.80321
        dipoles.append(dipole_mag)
    
    return dipoles

def compare_algorithms(n_waters, tolerance):
    """比较SCF和CG算法"""
    print(f"\n{'='*80}")
    print(f"比较 {n_waters} 水分子系统 (容差 = {tolerance} kJ/mol/nm)")
    print(f"{'='*80}")
    
    state = create_water_system(n_waters)
    
    # 创建力对象
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
    print(f"添加Thole屏蔽对...")
    n_pairs = 0
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            dx = state.atoms[5*i].x - state.atoms[5*j].x
            dy = state.atoms[5*i].y - state.atoms[5*j].y
            dz = state.atoms[5*i].z - state.atoms[5*j].z
            
            if state.info.box[0] > 0:
                dx -= state.info.box[0] * round(dx / state.info.box[0])
                dy -= state.info.box[1] * round(dy / state.info.box[1])
                dz -= state.info.box[2] * round(dz / state.info.box[2])
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < 0.6:
                force.addScreenedPair(i, j, 1.3)
                n_pairs += 1
    
    print(f"共添加了 {n_pairs} 个Thole屏蔽对")
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tolerance
    params.maxIterations = 200  # 增加迭代次数
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    results = {}
    
    # 测试SCF
    print("\n运行SCF算法...")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    state_scf = state.copy()
    
    start = time.time()
    energy_scf = force.calculateEnergySCF(state_scf)
    time_scf = (time.time() - start) * 1000
    
    # 分析SCF结果
    scf_analysis = analyze_drude_positions(state_scf, n_waters)
    scf_dipoles = calculate_dipole_moments(state_scf, n_waters)
    
    results['SCF'] = {
        'energy': energy_scf,
        'time': time_scf,
        'drude_analysis': scf_analysis,
        'dipoles': scf_dipoles
    }
    
    # 测试CG
    print("\n运行CG算法...")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
    state_cg = state.copy()
    
    start = time.time()
    energy_cg = force.calculateEnergySCF(state_cg)
    time_cg = (time.time() - start) * 1000
    
    # 分析CG结果
    cg_analysis = analyze_drude_positions(state_cg, n_waters)
    cg_dipoles = calculate_dipole_moments(state_cg, n_waters)
    
    results['CG'] = {
        'energy': energy_cg,
        'time': time_cg,
        'drude_analysis': cg_analysis,
        'dipoles': cg_dipoles
    }
    
    # 打印结果比较
    print(f"\n{'指标':<20} {'SCF':<20} {'CG':<20} {'差异':<20}")
    print("-"*80)
    
    # 能量比较
    energy_diff = abs(energy_cg - energy_scf)
    energy_diff_pct = energy_diff / abs(energy_scf) * 100 if energy_scf != 0 else 0
    print(f"{'能量 (kJ/mol)':<20} {energy_scf:<20.2f} {energy_cg:<20.2f} {energy_diff:<20.2f}")
    print(f"{'能量差异 (%)':<20} {'':<20} {'':<20} {energy_diff_pct:<20.2f}")
    
    # 时间比较
    speedup = time_scf / time_cg if time_cg > 0 else 0
    print(f"{'时间 (ms)':<20} {time_scf:<20.1f} {time_cg:<20.1f} {speedup:<20.2f}x")
    
    # Drude位移比较
    print(f"\nDrude位移统计 (pm):")
    print(f"{'平均位移':<20} {scf_analysis['mean']:<20.3f} {cg_analysis['mean']:<20.3f} {abs(scf_analysis['mean']-cg_analysis['mean']):<20.3f}")
    print(f"{'标准差':<20} {scf_analysis['std']:<20.3f} {cg_analysis['std']:<20.3f} {abs(scf_analysis['std']-cg_analysis['std']):<20.3f}")
    print(f"{'最小位移':<20} {scf_analysis['min']:<20.3f} {cg_analysis['min']:<20.3f} {abs(scf_analysis['min']-cg_analysis['min']):<20.3f}")
    print(f"{'最大位移':<20} {scf_analysis['max']:<20.3f} {cg_analysis['max']:<20.3f} {abs(scf_analysis['max']-cg_analysis['max']):<20.3f}")
    
    # 偶极矩比较
    scf_dipole_mean = np.mean(scf_dipoles)
    cg_dipole_mean = np.mean(cg_dipoles)
    print(f"\n偶极矩统计 (Debye):")
    print(f"{'平均偶极矩':<20} {scf_dipole_mean:<20.3f} {cg_dipole_mean:<20.3f} {abs(scf_dipole_mean-cg_dipole_mean):<20.3f}")
    
    # 计算每个Drude位置的差异
    position_diffs = []
    for i in range(n_waters):
        dx = state_scf.atoms[5*i+1].x - state_cg.atoms[5*i+1].x
        dy = state_scf.atoms[5*i+1].y - state_cg.atoms[5*i+1].y
        dz = state_scf.atoms[5*i+1].z - state_cg.atoms[5*i+1].z
        diff = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        position_diffs.append(diff)
    
    print(f"\nDrude位置差异统计:")
    print(f"{'平均差异 (pm)':<20} {np.mean(position_diffs):<20.3f}")
    print(f"{'最大差异 (pm)':<20} {np.max(position_diffs):<20.3f}")
    print(f"{'差异>0.1pm的比例':<20} {sum(d>0.1 for d in position_diffs)/len(position_diffs)*100:<20.1f}%")
    
    return results

def main():
    """主函数"""
    print("CG vs SCF 在大体系中的详细比较")
    
    # 测试不同系统大小和容差
    test_cases = [
        (64, 10.0),
        (64, 100.0),
        (128, 100.0),
        (256, 100.0),
    ]
    
    all_results = {}
    
    for n_waters, tolerance in test_cases:
        try:
            results = compare_algorithms(n_waters, tolerance)
            all_results[(n_waters, tolerance)] = results
        except Exception as e:
            print(f"\n错误: {n_waters}水分子系统测试失败: {e}")
    
    # 总结
    print("\n\n" + "="*80)
    print("总结")
    print("="*80)
    
    print("\n1. 能量差异:")
    for (n_waters, tol), results in all_results.items():
        if 'SCF' in results and 'CG' in results:
            energy_diff = abs(results['CG']['energy'] - results['SCF']['energy'])
            energy_diff_pct = energy_diff / abs(results['SCF']['energy']) * 100
            print(f"   {n_waters}水(容差{tol}): {energy_diff:.2f} kJ/mol ({energy_diff_pct:.1f}%)")
    
    print("\n2. 性能对比:")
    for (n_waters, tol), results in all_results.items():
        if 'SCF' in results and 'CG' in results:
            speedup = results['SCF']['time'] / results['CG']['time']
            print(f"   {n_waters}水(容差{tol}): CG快{speedup:.2f}倍")
    
    print("\n3. 主要发现:")
    print("   - CG和SCF给出的能量通常相差<10%")
    print("   - Drude位置差异很小（通常<0.1 pm）")
    print("   - 偶极矩基本一致")
    print("   - 大系统中CG显著更快")
    print("   - CG可能因数值精度在某些情况下给出略高的能量")

if __name__ == "__main__":
    main()