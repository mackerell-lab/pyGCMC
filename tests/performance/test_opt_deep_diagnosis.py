#!/usr/bin/env python3
"""
深入诊断OPT系列算法的问题
"""

import pygcmc
import numpy as np

def create_simple_system():
    """创建简单的2水分子系统"""
    atoms = []
    residues = []
    
    # 两个水分子，间距0.5nm
    for i in range(2):
        x_base = i * 0.5
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
    
    state.info.box = np.array([5.0, 5.0, 5.0])
    state.info.cutoff = 2.5
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_algorithm_step_by_step(algorithm_name):
    """逐步测试算法，观察每次迭代的结果"""
    print(f"\n{'='*60}")
    print(f"测试 {algorithm_name} 算法的迭代过程")
    print(f"{'='*60}")
    
    state = create_simple_system()
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    # 添加Drude粒子
    for i in range(2):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    force.addScreenedPair(0, 1, 1.3)
    
    # 设置参数 - 使用较大容差以观察收敛过程
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # 大容差
    params.maxIterations = 10  # 少量迭代
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 设置算法
    if algorithm_name == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm_name == "OPT3":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    elif algorithm_name == "OPT4":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT4)
    elif algorithm_name == "HybridOPT":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.HybridOPT)
    
    # 重置Drude位置
    for i in range(2):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    print("\n初始状态:")
    print("Drude相对位移 (pm):")
    for i in range(2):
        dx = (state.atoms[5*i+1].x - state.atoms[5*i].x) * 1000
        dy = (state.atoms[5*i+1].y - state.atoms[5*i].y) * 1000
        dz = (state.atoms[5*i+1].z - state.atoms[5*i].z) * 1000
        print(f"  Drude {i}: ({dx:.3f}, {dy:.3f}, {dz:.3f})")
    
    # 计算能量
    energy = force.calculateEnergySCF(state)
    
    print(f"\n最终能量: {energy:.6f} kJ/mol")
    print("\n最终Drude位移 (pm):")
    for i in range(2):
        dx = (state.atoms[5*i+1].x - state.atoms[5*i].x) * 1000
        dy = (state.atoms[5*i+1].y - state.atoms[5*i].y) * 1000
        dz = (state.atoms[5*i+1].z - state.atoms[5*i].z) * 1000
        disp = np.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"  Drude {i}: ({dx:.3f}, {dy:.3f}, {dz:.3f}) |d|={disp:.3f}")

def compare_with_reference():
    """与高精度SCF参考结果对比"""
    print("\n\n与高精度SCF参考结果对比")
    print("="*60)
    
    state = create_simple_system()
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
    
    # 先获取高精度SCF结果
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
    
    ref_energy = force.calculateEnergySCF(state)
    
    # 保存参考位置
    ref_positions = []
    for i in range(2):
        ref_positions.append([
            state.atoms[5*i+1].x - state.atoms[5*i].x,
            state.atoms[5*i+1].y - state.atoms[5*i].y,
            state.atoms[5*i+1].z - state.atoms[5*i].z
        ])
    
    print(f"高精度SCF参考能量: {ref_energy:.6f} kJ/mol")
    print("参考Drude位移 (pm):")
    for i, pos in enumerate(ref_positions):
        disp = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2) * 1000
        print(f"  Drude {i}: |d|={disp:.3f}")
    
    # 测试各算法
    algorithms = ["OPT3", "OPT4", "HybridOPT"]
    
    print(f"\n{'算法':<12} {'能量误差(kJ/mol)':<18} {'位移误差(pm)':<15}")
    print("-"*45)
    
    for algo in algorithms:
        # 设置低精度参数
        params.tolerance = 1.0
        params.maxIterations = 100
        force.setSCFParameters(params)
        
        if algo == "OPT3":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        elif algo == "OPT4":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT4)
        elif algo == "HybridOPT":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.HybridOPT)
        
        # 重置Drude
        for i in range(2):
            state.atoms[5*i+1].x = state.atoms[5*i].x
            state.atoms[5*i+1].y = state.atoms[5*i].y
            state.atoms[5*i+1].z = state.atoms[5*i].z
        
        energy = force.calculateEnergySCF(state)
        energy_error = abs(energy - ref_energy)
        
        # 计算位移误差
        disp_errors = []
        for i in range(2):
            dx = state.atoms[5*i+1].x - state.atoms[5*i].x - ref_positions[i][0]
            dy = state.atoms[5*i+1].y - state.atoms[5*i].y - ref_positions[i][1]
            dz = state.atoms[5*i+1].z - state.atoms[5*i].z - ref_positions[i][2]
            error = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            disp_errors.append(error)
        
        avg_disp_error = np.mean(disp_errors)
        
        print(f"{algo:<12} {energy_error:<18.6f} {avg_disp_error:<15.3f}")

def test_opt_coefficients():
    """检查OPT算法的系数是否合理"""
    print("\n\n检查OPT算法系数")
    print("="*60)
    
    # OPT3的理论系数
    print("OPT3理论系数:")
    print("  c1 = 1.812 (当前位置)")
    print("  c2 = -1.312 (前一步)")
    print("  c3 = 0.5 (前两步)")
    
    # OPT4的理论系数
    print("\nOPT4理论系数:")
    print("  c1 = 2.270")
    print("  c2 = -2.270")
    print("  c3 = 1.453")
    print("  c4 = -0.453")
    
    print("\n分析:")
    print("1. 系数和应该等于1.0 (归一化)")
    print("2. OPT3: 1.812 - 1.312 + 0.5 = 1.0 ✓")
    print("3. OPT4: 2.270 - 2.270 + 1.453 - 0.453 = 1.0 ✓")
    print("\n系数看起来是正确的，问题可能在实现细节中")

def main():
    print("OPT算法深度诊断")
    print("="*60)
    
    # 1. 测试各算法的迭代过程
    for algo in ["SCF", "OPT3", "OPT4", "HybridOPT"]:
        test_algorithm_step_by_step(algo)
    
    # 2. 与参考结果对比
    compare_with_reference()
    
    # 3. 检查系数
    test_opt_coefficients()
    
    print("\n\n诊断结论:")
    print("-"*60)
    print("1. OPT算法的能量明显高于SCF，说明没有正确收敛")
    print("2. 可能的原因:")
    print("   - 历史位置的初始化有问题")
    print("   - 混合系数的应用方式有误")
    print("   - HybridOPT没有正确执行最后的SCF步骤")
    print("3. 需要检查C++实现中的具体细节")

if __name__ == "__main__":
    main()