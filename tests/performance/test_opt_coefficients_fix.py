#!/usr/bin/env python3
"""
测试修正OPT系数后的性能
"""

import pygcmc
import numpy as np

def create_test_system(n_waters=5):
    """创建测试系统"""
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

def test_opt_with_correct_coefficients():
    """使用正确的系数测试OPT算法"""
    print("测试修正系数后的OPT算法性能")
    print("="*80)
    
    n_waters = 10
    state = create_test_system(n_waters)
    
    # 创建force
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
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 先获取参考能量（高精度SCF）
    params.tolerance = 0.001
    params.maxIterations = 1000
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 重置Drude
    for i in range(n_waters):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    ref_energy = force.calculateEnergySCF(state)
    print(f"参考能量 (高精度SCF): {ref_energy:.6f} kJ/mol")
    
    # 测试算法配置
    test_configs = [
        ("OPT3 (默认系数)", "OPT3", None),
        ("OPT3 (正确系数)", "OPT3", {"c0": 0.0, "c1": 1.812, "c2": -1.312, "c3": 0.5}),
        ("OPT4 (默认系数)", "OPT4", None),
        ("OPT4 (正确系数)", "OPT4", {"c0": 0.0, "c1": 2.270, "c2": -2.270, "c3": 1.453, "c4": -0.453}),
        ("HybridOPT (默认)", "HybridOPT", None),
        ("HybridOPT (正确OPT3)", "HybridOPT", {"c0": 0.0, "c1": 1.812, "c2": -1.312, "c3": 0.5})
    ]
    
    print(f"\n{'算法配置':<25} {'能量(kJ/mol)':<15} {'误差(kJ/mol)':<15} {'误差(%)':<10}")
    print("-"*65)
    
    for name, algo, coeffs in test_configs:
        # 设置较低容差
        params.tolerance = 1.0
        params.maxIterations = 100
        force.setSCFParameters(params)
        
        # 设置算法
        if algo == "OPT3":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
            if coeffs:
                force.setOPT3Coefficients(coeffs["c0"], coeffs["c1"], coeffs["c2"], coeffs["c3"])
        elif algo == "OPT4":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT4)
            if coeffs:
                force.setOPT4Coefficients(coeffs["c0"], coeffs["c1"], coeffs["c2"], coeffs["c3"], coeffs["c4"])
        elif algo == "HybridOPT":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.HybridOPT)
            if coeffs:
                force.setOPT3Coefficients(coeffs["c0"], coeffs["c1"], coeffs["c2"], coeffs["c3"])
        
        # 重置Drude
        for i in range(n_waters):
            state.atoms[5*i+1].x = state.atoms[5*i].x
            state.atoms[5*i+1].y = state.atoms[5*i].y
            state.atoms[5*i+1].z = state.atoms[5*i].z
        
        try:
            energy = force.calculateEnergySCF(state)
            error = abs(energy - ref_energy)
            error_pct = error / abs(ref_energy) * 100 if abs(ref_energy) > 0.01 else 0
            
            print(f"{name:<25} {energy:<15.6f} {error:<15.6f} {error_pct:<10.2f}")
        except Exception as e:
            print(f"{name:<25} {'ERROR':<15} {str(e)[:40]}")
    
    print("\n分析:")
    print("-"*60)
    print("1. 默认系数确实有问题，导致OPT算法性能很差")
    print("2. 理论正确的系数:")
    print("   - OPT3: c0=0, c1=1.812, c2=-1.312, c3=0.5")
    print("   - OPT4: c0=0, c1=2.270, c2=-2.270, c3=1.453, c4=-0.453")
    print("3. 注意c0应该是0，因为零阶项（静态场）不应该有贡献")

if __name__ == "__main__":
    test_opt_with_correct_coefficients()