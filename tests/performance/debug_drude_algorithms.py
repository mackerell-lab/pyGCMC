#!/usr/bin/env python3
"""
调试为什么其他Drude算法不工作
"""

import numpy as np
import pygcmc

def test_each_algorithm():
    """
    逐个测试每种算法
    """
    print("调试Drude算法")
    print("="*70)
    
    # 创建最简单的2水系统
    state = pygcmc.MCState()
    
    box_size = 1.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    # 两个水分子
    positions = [
        # 水1
        [0.3, 0.3, 0.3],      # O
        [0.3, 0.3, 0.3],      # D
        [0.396, 0.3, 0.3],    # H1
        [0.252, 0.377, 0.3],  # H2
        [0.3, 0.3, 0.3],      # M
        # 水2
        [0.7, 0.7, 0.7],      # O
        [0.7, 0.7, 0.7],      # D
        [0.796, 0.7, 0.7],    # H1
        [0.652, 0.777, 0.7],  # H2
        [0.7, 0.7, 0.7]       # M
    ]
    
    atoms = []
    residues = []
    
    for i in range(10):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = i % 5 if i % 5 < 4 else 3
        atoms.append(atom)
    
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 10
    state.activeResidueCount = 2
    
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print("测试系统：2个水分子")
    
    # 测试每种算法
    algorithms = [
        (pygcmc.DrudeAlgorithm.SCF, "SCF", "calculateEnergySCF"),
        (pygcmc.DrudeAlgorithm.ConjugateGradient, "ConjugateGradient", "calculateEnergyCG"),
        (pygcmc.DrudeAlgorithm.OPT3, "OPT3", "calculateEnergyOPT"),
        (pygcmc.DrudeAlgorithm.OPT4, "OPT4", "calculateEnergyOPT"),
        (pygcmc.DrudeAlgorithm.SmartOPT3, "SmartOPT3", "calculateEnergyOPT"),
        (pygcmc.DrudeAlgorithm.AdaptiveOPT, "AdaptiveOPT", "calculateEnergyOPT"),
        (pygcmc.DrudeAlgorithm.HybridOPT, "HybridOPT", "calculateEnergyOPT"),
        (pygcmc.DrudeAlgorithm.FBP, "FBP", "calculateEnergyFBP")
    ]
    
    print("\n测试结果：")
    print(f"{'算法':>20} {'状态':>10} {'能量(kJ/mol)':>15} {'错误信息':>40}")
    print("-"*90)
    
    for algo, name, method_name in algorithms:
        # 创建新的DrudeForce
        force = pygcmc.DrudeForce()
        
        # 添加Drude粒子
        for i in range(2):
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
        
        # 添加Thole对
        force.addScreenedPair(0, 1, 1.3)
        
        # 设置算法
        force.setAlgorithm(algo)
        
        # 对于SCF设置参数
        if algo == pygcmc.DrudeAlgorithm.SCF:
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 10.0
            params.maxIterations = 100
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
        
        # 创建测试状态
        state_test = state.copy()
        
        # 测试
        try:
            # 获取正确的方法
            if hasattr(force, method_name):
                method = getattr(force, method_name)
                energy = method(state_test)
                
                # 计算平均位移
                displacements = []
                for i in range(2):
                    o_idx = i * 5
                    d_idx = i * 5 + 1
                    
                    dx = state_test.atoms[d_idx].x - state_test.atoms[o_idx].x
                    dy = state_test.atoms[d_idx].y - state_test.atoms[o_idx].y
                    dz = state_test.atoms[d_idx].z - state_test.atoms[o_idx].z
                    
                    disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                    displacements.append(disp)
                
                avg_disp = np.mean(displacements)
                
                print(f"{name:>20} {'成功':>10} {energy:>15.4f} {f'位移={avg_disp:.2f}pm':>40}")
            else:
                print(f"{name:>20} {'失败':>10} {'-':>15} {f'方法{method_name}不存在':>40}")
                
        except Exception as e:
            error_msg = str(e)[:40]
            print(f"{name:>20} {'异常':>10} {'-':>15} {error_msg:>40}")
    
    # 列出force的所有方法
    print("\n\nDrudeForce可用的方法：")
    methods = [m for m in dir(force) if 'calculate' in m.lower() and not m.startswith('_')]
    for method in sorted(methods):
        print(f"  - {method}")

def test_algorithm_requirements():
    """
    测试不同算法的要求
    """
    print("\n\n测试算法要求")
    print("="*70)
    
    # 检查是否需要特殊设置
    print("\n可能的原因：")
    print("1. 某些算法可能需要额外的参数设置")
    print("2. 某些算法可能需要特定的系统配置")
    print("3. 某些算法可能还未完全实现")
    print("4. 方法名称可能不同")
    
    # 测试是否有统一的接口
    print("\n测试统一接口...")
    
    force = pygcmc.DrudeForce()
    force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=-1.71636,
        polarizability=0.0009782237,
        aniso12=1.0,
        aniso34=1.0
    )
    
    # 检查是否有通用的计算方法
    if hasattr(force, 'calculateEnergy'):
        print("  ✓ 找到calculateEnergy方法")
    else:
        print("  ✗ 没有通用的calculateEnergy方法")

def main():
    """
    主函数
    """
    test_each_algorithm()
    test_algorithm_requirements()
    
    print("\n\n结论：")
    print("="*70)
    print("1. 不同算法可能需要不同的方法调用")
    print("2. SCF是最稳定和成熟的实现")
    print("3. 其他算法可能还在开发中")

if __name__ == "__main__":
    main()