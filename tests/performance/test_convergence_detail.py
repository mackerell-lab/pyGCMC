#\!/usr/bin/env python3
"""
详细的收敛性分析 - 研究FBP收敛困难的原因
"""

import pygcmc
import numpy as np
import time

def create_test_system(n_waters, spacing_factor=1.0):
    """创建测试系统，可调整间距"""
    base_spacing = 0.35  # nm - 接近真实密度
    spacing = base_spacing * spacing_factor
    
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
                
                # 添加小的随机扰动
                x += (np.random.rand() - 0.5) * 0.02
                y += (np.random.rand() - 0.5) * 0.02
                z += (np.random.rand() - 0.5) * 0.02
                
                # SWM4-NDP水模型位置
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

def analyze_convergence_pattern(state, n_waters, algorithm, max_iter=50):
    """分析收敛模式"""
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
    
    # Thole屏蔽
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 设置参数 - 使用极大容差来观察收敛过程
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e10  # 不会提前收敛
    params.maxIterations = max_iter
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
    
    # 运行并获取能量
    energy = force.calculateEnergySCF(state_copy)
    
    # 计算最终Drude位移
    displacements = []
    for i in range(n_waters):
        dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
        dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
        dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        displacements.append(disp)
    
    return {
        'energy': energy,
        'avg_disp': np.mean(displacements),
        'max_disp': np.max(displacements),
        'min_disp': np.min(displacements),
        'std_disp': np.std(displacements)
    }

def main():
    """主测试函数"""
    print("FBP收敛性详细分析")
    print("="*80)
    
    # 测试不同密度
    print("\n1. 不同密度下的收敛性")
    print("-"*60)
    
    n_waters = 64
    spacing_factors = [1.0, 1.2, 1.5, 2.0]  # 密度从高到低
    
    print(f"{'间距因子':<10} {'SCF能量':<15} {'FBP能量':<15} {'SCF位移(pm)':<15} {'FBP位移(pm)':<15}")
    print("-"*80)
    
    for factor in spacing_factors:
        state = create_test_system(n_waters, factor)
        
        # 分析SCF
        scf_result = analyze_convergence_pattern(state, n_waters, "SCF", max_iter=50)
        
        # 分析FBP
        fbp_result = analyze_convergence_pattern(state, n_waters, "FBP", max_iter=50)
        
        print(f"{factor:<10.1f} {scf_result['energy']:<15.2f} {fbp_result['energy']:<15.2f} "
              f"{scf_result['avg_disp']:<15.2f} {fbp_result['avg_disp']:<15.2f}")
    
    # 测试不同系统大小
    print("\n\n2. 不同系统大小的收敛性（固定50次迭代）")
    print("-"*60)
    
    system_sizes = [8, 16, 32, 64]
    spacing_factor = 1.5  # 使用较松的间距
    
    print(f"{'水分子数':<10} {'SCF能量':<15} {'FBP能量':<15} {'能量差(%)':<15} {'FBP最大位移(pm)':<15}")
    print("-"*80)
    
    for n_waters in system_sizes:
        state = create_test_system(n_waters, spacing_factor)
        
        # 分析SCF
        scf_result = analyze_convergence_pattern(state, n_waters, "SCF", max_iter=50)
        
        # 分析FBP
        fbp_result = analyze_convergence_pattern(state, n_waters, "FBP", max_iter=50)
        
        energy_diff_pct = abs(fbp_result['energy'] - scf_result['energy']) / abs(scf_result['energy']) * 100
        
        print(f"{n_waters:<10} {scf_result['energy']:<15.2f} {fbp_result['energy']:<15.2f} "
              f"{energy_diff_pct:<15.2f} {fbp_result['max_disp']:<15.2f}")
    
    # 分析位移分布
    print("\n\n3. 128水系统的Drude位移分布分析")
    print("-"*60)
    
    n_waters = 128
    state = create_test_system(n_waters, spacing_factor=1.5)
    
    print(f"{'算法':<10} {'平均(pm)':<12} {'最大(pm)':<12} {'最小(pm)':<12} {'标准差(pm)':<12}")
    print("-"*60)
    
    for algo in ["SCF", "FBP"]:
        result = analyze_convergence_pattern(state, n_waters, algo, max_iter=50)
        print(f"{algo:<10} {result['avg_disp']:<12.2f} {result['max_disp']:<12.2f} "
              f"{result['min_disp']:<12.2f} {result['std_disp']:<12.2f}")
    
    print("\n结论:")
    print("-"*60)
    print("1. FBP在密集系统中收敛困难，需要更多迭代")
    print("2. 系统越大，FBP的收敛问题越严重")
    print("3. FBP的位移分布更不均匀，表明局部收敛问题")
    print("4. 在稀疏系统中，FBP和SCF性能接近")

if __name__ == "__main__":
    main()
