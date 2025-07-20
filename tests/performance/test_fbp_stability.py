#!/usr/bin/env python3
"""
测试FBP在不同配置下的稳定性和准确性
"""

import pygcmc
import numpy as np

def create_water_cluster(n_waters, spacing):
    """创建水分子团簇"""
    atoms = []
    residues = []
    
    # 将水分子排列在立方体顶点
    positions = []
    if n_waters == 1:
        positions = [(0, 0, 0)]
    elif n_waters == 2:
        positions = [(0, 0, 0), (spacing, 0, 0)]
    elif n_waters == 4:
        positions = [
            (0, 0, 0), (spacing, 0, 0),
            (0, spacing, 0), (spacing, spacing, 0)
        ]
    elif n_waters == 8:
        positions = [
            (0, 0, 0), (spacing, 0, 0),
            (0, spacing, 0), (spacing, spacing, 0),
            (0, 0, spacing), (spacing, 0, spacing),
            (0, spacing, spacing), (spacing, spacing, spacing)
        ]
    
    for i, (x_base, y_base, z_base) in enumerate(positions[:n_waters]):
        # 水分子原子
        water_atoms = [
            (x_base, y_base, z_base, 1.71636, 0),   # O
            (x_base, y_base, z_base, -1.71636, 1),  # D
            (x_base + 0.09572, y_base, z_base, 0.55733, 2),  # H1
            (x_base - 0.04786, y_base + 0.08288, z_base, 0.55733, 2),  # H2
            (x_base, y_base - 0.024034, z_base, -1.11466, 3)  # M
        ]
        
        for x, y, z, charge, typ in water_atoms:
            atom = pygcmc.MCAtom()
            atom.x = x
            atom.y = y
            atom.z = z
            atom.charge = charge
            atom.type = typ
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    # 盒子稍大于系统
    box_size = max(10.0, spacing * 2 + 2.0)
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(4.5, box_size/2 - 0.1)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_algorithm(state, n_waters, algorithm, tolerance=0.1, verbose=False):
    """测试算法并返回详细结果"""
    force = pygcmc.DrudeForce()
    
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    # 添加Drude粒子
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 添加屏蔽对
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    if algorithm == "SCF":
        params.tolerance = 0.001
        params.maxIterations = 500
    else:
        params.tolerance = tolerance
        params.maxIterations = 50
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm == "FBP":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    elif algorithm == "OPT3":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    
    # 重置Drude位置
    for i in range(n_waters):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    # 计算能量
    energy = force.calculateEnergySCF(state)
    
    # 收集Drude位移
    displacements = []
    for i in range(n_waters):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        dy = state.atoms[5*i+1].y - state.atoms[5*i].y
        dz = state.atoms[5*i+1].z - state.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz)
        displacements.append(disp)
    
    avg_disp = np.mean(displacements) * 1000  # nm to pm
    max_disp = np.max(displacements) * 1000
    
    if verbose:
        print(f"  平均位移: {avg_disp:.3f} pm")
        print(f"  最大位移: {max_disp:.3f} pm")
        print(f"  位移分布: ", end="")
        for d in displacements[:5]:
            print(f"{d*1000:.3f} ", end="")
        if len(displacements) > 5:
            print("...")
        else:
            print()
    
    return energy, avg_disp, max_disp, displacements

def main():
    """主测试函数"""
    print("FBP稳定性和准确性测试")
    print("="*80)
    
    # 测试不同的水分子间距
    spacings = [0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
    n_waters_list = [1, 2, 4, 8]
    
    print("\n1. 不同间距下的能量比较")
    print("-"*80)
    
    for n_waters in [2, 4]:
        print(f"\n{n_waters}个水分子:")
        print(f"{'Spacing (nm)':<15} {'SCF Energy':<15} {'FBP Energy':<15} {'Error (%)':<15} {'FBP Disp (pm)':<15}")
        print("-"*75)
        
        for spacing in spacings:
            if n_waters > 2 and spacing < 0.4:
                continue  # 跳过太拥挤的配置
                
            state = create_water_cluster(n_waters, spacing)
            
            # SCF参考
            energy_scf, _, _, _ = test_algorithm(state, n_waters, "SCF")
            
            # FBP测试
            energy_fbp, avg_disp_fbp, _, _ = test_algorithm(state, n_waters, "FBP")
            
            # 计算误差
            if abs(energy_scf) > 0.01:
                error = abs(energy_fbp - energy_scf) / abs(energy_scf) * 100
            else:
                error = abs(energy_fbp - energy_scf) * 100
            
            print(f"{spacing:<15.2f} {energy_scf:<15.6f} {energy_fbp:<15.6f} {error:<15.2f} {avg_disp_fbp:<15.3f}")
    
    # 测试FBP的收敛容差影响
    print("\n\n2. FBP容差对精度的影响 (4个水分子，间距0.5nm)")
    print("-"*80)
    
    state = create_water_cluster(4, 0.5)
    energy_scf, _, _, _ = test_algorithm(state, 4, "SCF")
    
    print(f"SCF参考能量: {energy_scf:.6f} kJ/mol")
    print(f"\n{'Tolerance':<15} {'FBP Energy':<15} {'Error':<15} {'Error (%)':<15}")
    print("-"*60)
    
    tolerances = [10.0, 5.0, 2.0, 1.0, 0.5, 0.1, 0.05, 0.01]
    for tol in tolerances:
        energy_fbp, _, _, _ = test_algorithm(state, 4, "FBP", tolerance=tol)
        error = abs(energy_fbp - energy_scf)
        error_pct = error / abs(energy_scf) * 100 if abs(energy_scf) > 0.01 else error * 100
        print(f"{tol:<15.2f} {energy_fbp:<15.6f} {error:<15.6f} {error_pct:<15.2f}")
    
    # 分析FBP失败的情况
    print("\n\n3. 详细分析FBP误差较大的情况")
    print("-"*80)
    
    # 选择一个误差较大的配置
    state = create_water_cluster(8, 0.3)
    n_waters = 8
    
    print(f"\n密集系统 ({n_waters}个水分子，间距0.3nm):")
    
    # SCF
    print("\nSCF结果:")
    energy_scf, avg_disp_scf, max_disp_scf, disp_scf = test_algorithm(
        state, n_waters, "SCF", verbose=True
    )
    print(f"  能量: {energy_scf:.6f} kJ/mol")
    
    # FBP
    print("\nFBP结果:")
    energy_fbp, avg_disp_fbp, max_disp_fbp, disp_fbp = test_algorithm(
        state, n_waters, "FBP", verbose=True
    )
    print(f"  能量: {energy_fbp:.6f} kJ/mol")
    print(f"  能量误差: {abs(energy_fbp - energy_scf):.6f} kJ/mol ({abs(energy_fbp - energy_scf)/abs(energy_scf)*100:.1f}%)")
    
    # OPT3比较
    print("\nOPT3结果:")
    energy_opt3, avg_disp_opt3, max_disp_opt3, _ = test_algorithm(
        state, n_waters, "OPT3", verbose=True
    )
    print(f"  能量: {energy_opt3:.6f} kJ/mol")
    print(f"  能量误差: {abs(energy_opt3 - energy_scf):.6f} kJ/mol ({abs(energy_opt3 - energy_scf)/abs(energy_scf)*100:.1f}%)")
    
    # 总结
    print("\n\n总结:")
    print("="*60)
    print("1. FBP在水分子间距较大时精度较好")
    print("2. 在密集系统中FBP误差增大")
    print("3. FBP的误差可能来自于多体效应的处理")
    print("4. 力平衡不等于能量最小，可能收敛到了局部极小值")

if __name__ == "__main__":
    main()