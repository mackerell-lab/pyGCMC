#!/usr/bin/env python3
"""
测试不同算法优化后的能量是否一致
"""

import pygcmc
import numpy as np

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

def test_energy_consistency():
    """测试不同算法的能量一致性"""
    print("能量一致性测试")
    print("="*80)
    print("\n说明：如果所有算法都正确，它们应该收敛到相同的最小能量")
    print("-"*80)
    
    # 测试不同大小的系统
    n_waters_list = [2, 5, 10]
    
    for n_waters in n_waters_list:
        print(f"\n\n{n_waters}水分子系统")
        print("="*60)
        
        state = create_water_system(n_waters)
        
        # 测试配置
        test_configs = [
            ("SCF (高精度)", "SCF", 0.001, None),
            ("SCF (中精度)", "SCF", 0.1, None),
            ("SCF (低精度)", "SCF", 1.0, None),
            ("FBP", "FBP", 1.0, None),
            ("OPT3 (默认系数)", "OPT3", 1.0, None),
            ("OPT3 (正确系数)", "OPT3", 1.0, {"c0": 0.0, "c1": 1.812, "c2": -1.312, "c3": 0.5}),
            ("SmartOPT3", "SmartOPT3", 1.0, None)
        ]
        
        results = []
        
        for name, algo, tol, coeffs in test_configs:
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
            
            # 设置参数
            params = pygcmc.DrudeSCFParams()
            params.tolerance = tol
            params.maxIterations = 1000
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            
            # 设置算法
            if algo == "SCF":
                force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
            elif algo == "FBP":
                force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
            elif algo == "OPT3":
                force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
                if coeffs:
                    force.setOPT3Coefficients(coeffs["c0"], coeffs["c1"], coeffs["c2"], coeffs["c3"])
            elif algo == "SmartOPT3":
                force.setAlgorithm(pygcmc.DrudeAlgorithm.SmartOPT3)
            
            # 重置Drude位置
            state_copy = state.copy()
            for i in range(n_waters):
                state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
                state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
                state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
            
            try:
                energy = force.calculateEnergySCF(state_copy)
                
                # 计算Drude位移
                displacements = []
                for i in range(n_waters):
                    dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
                    dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
                    dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
                    disp = np.sqrt(dx*dx + dy*dy + dz*dz)
                    displacements.append(disp)
                
                avg_disp = np.mean(displacements) * 1000  # nm to pm
                max_disp = np.max(displacements) * 1000
                
                results.append({
                    'name': name,
                    'energy': energy,
                    'avg_disp': avg_disp,
                    'max_disp': max_disp
                })
                
            except Exception as e:
                results.append({
                    'name': name,
                    'energy': float('nan'),
                    'avg_disp': float('nan'),
                    'max_disp': float('nan')
                })
        
        # 找到最低能量作为参考
        valid_energies = [r['energy'] for r in results if not np.isnan(r['energy'])]
        if valid_energies:
            min_energy = min(valid_energies)
            
            # 打印结果
            print(f"\n{'算法':<25} {'能量(kJ/mol)':<15} {'与最低能量差':<15} {'平均位移(pm)':<12} {'最大位移(pm)':<12}")
            print("-"*90)
            
            for r in results:
                if not np.isnan(r['energy']):
                    diff = r['energy'] - min_energy
                    print(f"{r['name']:<25} {r['energy']:<15.6f} {diff:<15.6f} {r['avg_disp']:<12.3f} {r['max_disp']:<12.3f}")
                else:
                    print(f"{r['name']:<25} {'失败':<15} {'N/A':<15} {'N/A':<12} {'N/A':<12}")
            
            # 分析能量分布
            print(f"\n能量分析:")
            print(f"  最低能量: {min_energy:.6f} kJ/mol")
            print(f"  最高能量: {max(valid_energies):.6f} kJ/mol")
            print(f"  能量范围: {max(valid_energies) - min_energy:.6f} kJ/mol")
            
            # 检查是否收敛到相同能量
            energy_std = np.std([r['energy'] for r in results if not np.isnan(r['energy'])])
            print(f"  能量标准差: {energy_std:.6f} kJ/mol")
            
            if energy_std < 0.01:
                print("  ✓ 所有算法收敛到相同能量（误差<0.01 kJ/mol）")
            else:
                print("  ✗ 算法未收敛到相同能量！")

def analyze_energy_components():
    """分析能量组成"""
    print("\n\n能量组成分析")
    print("="*80)
    
    # 创建简单的2水分子系统
    state = create_water_system(2)
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    for i in range(2):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    force.addScreenedPair(0, 1, 1.3)
    
    # 高精度SCF
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.001
    params.maxIterations = 1000
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 重置Drude
    for i in range(2):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    total_energy = force.calculateEnergySCF(state)
    
    print(f"总能量: {total_energy:.6f} kJ/mol")
    
    # 分析各项贡献
    print("\n能量组成（估算）:")
    print("1. Drude谐振子能量 (k/2 * r²)")
    print("2. 库仑相互作用能量")
    print("3. Thole屏蔽修正")
    
    # 计算Drude谐振子能量
    spring_energy = 0.0
    for i in range(2):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        dy = state.atoms[5*i+1].y - state.atoms[5*i].y
        dz = state.atoms[5*i+1].z - state.atoms[5*i].z
        r2 = dx*dx + dy*dy + dz*dz
        spring_energy += 0.5 * k_spring * r2
    
    print(f"\nDrude弹簧能量: {spring_energy:.6f} kJ/mol")
    print(f"其他相互作用: {total_energy - spring_energy:.6f} kJ/mol")
    
    # 为什么能量是正的？
    print("\n\n为什么能量是正的？")
    print("-"*60)
    print("1. 这是相对于'无限远分离的中性原子'的能量")
    print("2. 正能量表示系统处于排斥状态")
    print("3. 在0.5nm间距下，两个水分子确实是相互排斥的")
    print("4. 如果增加间距，能量会降低并可能变为负值")

if __name__ == "__main__":
    test_energy_consistency()
    analyze_energy_components()