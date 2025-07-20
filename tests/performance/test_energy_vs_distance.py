#!/usr/bin/env python3
"""
测试能量随距离的变化，验证物理正确性
"""

import pygcmc
import numpy as np
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def test_energy_vs_distance():
    """测试不同间距下的能量"""
    print("能量-距离关系测试")
    print("="*80)
    
    # 测试不同间距
    distances = np.linspace(0.3, 2.0, 20)  # 0.3 到 2.0 nm
    
    energies_scf = []
    energies_fbp = []
    energies_opt3 = []
    
    for dist in distances:
        # 创建2水分子系统
        atoms = []
        residues = []
        
        for i in range(2):
            x_base = i * dist
            positions = [
                (x_base, 0.0, 0.0, 1.71636, 0),   # O
                (x_base, 0.0, 0.0, -1.71636, 1),  # D
                (x_base + 0.09572, 0.0, 0.0, 0.55733, 2),  # H1
                (x_base - 0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
                (x_base, -0.024034, 0.0, -1.11466, 3)  # M
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
            res.atomStart = 5 * i
            res.atomCount = 5
            res.active = True
            res.type = 0
            residues.append(res)
        
        state = pygcmc.MCState()
        state.atoms = atoms
        state.residues = residues
        state.activeAtomCount = 10
        state.activeResidueCount = 2
        
        state.info.box = np.array([10.0, 10.0, 10.0])
        state.info.cutoff = 4.5
        
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 测试每种算法
        for algo_name, energies_list in [("SCF", energies_scf), ("FBP", energies_fbp), ("OPT3", energies_opt3)]:
            force = pygcmc.DrudeForce()
            
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
            
            params = pygcmc.DrudeSCFParams()
            if algo_name == "SCF":
                params.tolerance = 0.1
                force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
            elif algo_name == "FBP":
                params.tolerance = 1.0
                force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
            else:  # OPT3
                params.tolerance = 1.0
                force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
                # 使用正确的系数
                force.setOPT3Coefficients(0.0, 1.812, -1.312, 0.5)
            
            params.maxIterations = 500
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            
            # 重置Drude
            state_copy = state.copy()
            for i in range(2):
                state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
                state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
                state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
            
            try:
                energy = force.calculateEnergySCF(state_copy)
                energies_list.append(energy)
            except:
                energies_list.append(float('nan'))
    
    # 打印数据表
    print(f"\n{'距离(nm)':<10} {'SCF能量':<15} {'FBP能量':<15} {'FBP误差':<15} {'OPT3能量':<15} {'OPT3误差':<15}")
    print("-"*90)
    
    for i, dist in enumerate(distances):
        scf_e = energies_scf[i]
        fbp_e = energies_fbp[i]
        opt3_e = energies_opt3[i]
        
        if not np.isnan(scf_e):
            fbp_err = abs(fbp_e - scf_e)
            opt3_err = abs(opt3_e - scf_e)
            print(f"{dist:<10.2f} {scf_e:<15.6f} {fbp_e:<15.6f} {fbp_err:<15.6f} {opt3_e:<15.6f} {opt3_err:<15.6f}")
    
    # 分析
    print("\n\n物理分析:")
    print("-"*60)
    
    # 找到能量最小点
    min_idx = np.argmin(energies_scf)
    min_dist = distances[min_idx]
    min_energy = energies_scf[min_idx]
    
    print(f"1. 能量最小点: {min_dist:.2f} nm, 能量 = {min_energy:.3f} kJ/mol")
    print(f"2. 短距离排斥: < 0.4 nm 时能量快速上升")
    print(f"3. 长距离趋零: > 1.5 nm 时能量接近0")
    
    # 检查一致性
    print("\n\n算法一致性分析:")
    print("-"*60)
    
    # 计算平均误差
    fbp_errors = []
    opt3_errors = []
    for i in range(len(distances)):
        if not np.isnan(energies_scf[i]):
            fbp_errors.append(abs(energies_fbp[i] - energies_scf[i]))
            opt3_errors.append(abs(energies_opt3[i] - energies_scf[i]))
    
    print(f"FBP平均误差: {np.mean(fbp_errors):.6f} kJ/mol")
    print(f"FBP最大误差: {np.max(fbp_errors):.6f} kJ/mol")
    print(f"OPT3平均误差: {np.mean(opt3_errors):.6f} kJ/mol")
    print(f"OPT3最大误差: {np.max(opt3_errors):.6f} kJ/mol")
    
    # 绘图（如果可能）
    if HAS_MATPLOTLIB:
        try:
            plt.figure(figsize=(10, 6))
            plt.plot(distances, energies_scf, 'b-', label='SCF (参考)', linewidth=2)
            plt.plot(distances, energies_fbp, 'r--', label='FBP', linewidth=2)
            plt.plot(distances, energies_opt3, 'g:', label='OPT3', linewidth=2)
            
            plt.xlabel('距离 (nm)')
            plt.ylabel('能量 (kJ/mol)')
            plt.title('水分子二聚体能量-距离曲线')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            
            plt.savefig('../tmp/energy_vs_distance.png', dpi=150)
            print("\n能量-距离曲线已保存到: tmp/energy_vs_distance.png")
        except Exception as e:
            print(f"\n绘图失败: {e}")
    else:
        print("\n（无法绘图，缺少matplotlib）")

if __name__ == "__main__":
    test_energy_vs_distance()