#!/usr/bin/env python3
"""
快速比较Drude算法性能
"""

import numpy as np
import time
import pickle
import os
import pygcmc

def quick_benchmark():
    """
    快速性能测试
    """
    print("Drude算法快速基准测试")
    print("="*70)
    
    # 只测试256水系统
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    if not os.path.exists(filename):
        print("文件不存在")
        return
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = 256
    positions = data['positions']
    box_length = data['box_length']
    charges = data['charges']
    
    print(f"测试系统：{n_waters}水分子")
    print(f"盒子大小：{box_length:.3f} nm")
    
    # 创建状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters * 5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = atom_types[i % 5]
        atoms.append(atom)
    
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(0.9, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 测试的算法
    algorithms = [
        (pygcmc.DrudeAlgorithm.SCF, "SCF", {'tolerance': 100.0, 'maxIterations': 50}),
        (pygcmc.DrudeAlgorithm.ConjugateGradient, "CG", None),
        (pygcmc.DrudeAlgorithm.OPT3, "OPT3", None),
        (pygcmc.DrudeAlgorithm.SmartOPT3, "SmartOPT3", None),
        (pygcmc.DrudeAlgorithm.HybridOPT, "HybridOPT", None),
        (pygcmc.DrudeAlgorithm.FBP, "FBP", None)
    ]
    
    print("\n测试结果：")
    print(f"{'算法':>15} {'时间(ms)':>12} {'能量/水(kJ/mol)':>18} {'平均位移(pm)':>15} {'状态':>10}")
    print("-"*75)
    
    best_time = float('inf')
    best_algo = None
    
    for algo, algo_name, params in algorithms:
        # 创建DrudeForce
        force = pygcmc.DrudeForce()
        
        # 添加Drude粒子
        for i in range(n_waters):
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
        
        # 添加少量Thole对（加速测试）
        n_thole = 0
        for i in range(min(50, n_waters)):
            for j in range(i+1, min(i+10, n_waters)):
                force.addScreenedPair(i, j, 1.3)
                n_thole += 1
        
        # 设置算法
        force.setAlgorithm(algo)
        
        # 设置参数
        if params and algo == pygcmc.DrudeAlgorithm.SCF:
            scf_params = pygcmc.DrudeSCFParams()
            scf_params.tolerance = params['tolerance']
            scf_params.maxIterations = params['maxIterations']
            scf_params.dampingFactor = 0.5
            scf_params.maxDrudeDistance = 0.02
            force.setSCFParameters(scf_params)
        
        # 准备测试状态
        state_test = state.copy()
        for i in range(n_waters):
            o_idx = i * 5
            d_idx = i * 5 + 1
            state_test.atoms[d_idx].x = state_test.atoms[o_idx].x
            state_test.atoms[d_idx].y = state_test.atoms[o_idx].y
            state_test.atoms[d_idx].z = state_test.atoms[o_idx].z
        
        # 测试
        try:
            start_time = time.time()
            
            if algo == pygcmc.DrudeAlgorithm.SCF:
                energy = force.calculateEnergySCF(state_test)
            elif algo == pygcmc.DrudeAlgorithm.ConjugateGradient:
                energy = force.calculateEnergyCG(state_test)
            else:
                energy = force.calculateEnergyOPT(state_test)
            
            end_time = time.time()
            elapsed_time = (end_time - start_time) * 1000  # ms
            
            # 计算平均位移
            displacements = []
            for i in range(min(50, n_waters)):  # 只检查前50个
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = state_test.atoms[d_idx].x - state_test.atoms[o_idx].x
                dy = state_test.atoms[d_idx].y - state_test.atoms[o_idx].y
                dz = state_test.atoms[d_idx].z - state_test.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                displacements.append(disp)
            
            avg_disp = np.mean(displacements)
            
            print(f"{algo_name:>15} {elapsed_time:>12.2f} {energy/n_waters:>18.2f} {avg_disp:>15.2f} {'成功':>10}")
            
            if elapsed_time < best_time:
                best_time = elapsed_time
                best_algo = algo_name
            
        except Exception as e:
            print(f"{algo_name:>15} {'-':>12} {'-':>18} {'-':>15} {'失败':>10}")
            # print(f"  错误: {str(e)[:50]}")
    
    print("\n总结：")
    print(f"最快算法：{best_algo} ({best_time:.2f} ms)")
    
    # 性能分析
    print("\n\n详细分析：")
    print("="*70)
    
    print("\n1. 算法特点：")
    print("   - SCF：迭代求解，精度高但可能慢")
    print("   - CG：共轭梯度法，适合大系统")
    print("   - OPT3：优化的直接求解，速度快")
    print("   - SmartOPT3：自适应版本")
    print("   - HybridOPT：混合方法，可能结合多种算法")
    print("   - FBP：力平衡方法，简单快速")
    
    print("\n2. 选择建议：")
    print("   - 精度要求高：使用SCF")
    print("   - 速度要求高：使用OPT3或FBP")
    print("   - 大系统：使用CG或FBP")
    print("   - 一般用途：SmartOPT3或HybridOPT")

def test_scaling():
    """
    测试算法随系统大小的扩展性
    """
    print("\n\n扩展性测试")
    print("="*70)
    
    # 创建不同大小的测试系统
    system_sizes = [50, 100, 200]
    
    print(f"\n{'系统大小':>10} {'SCF(ms)':>10} {'OPT3(ms)':>10} {'FBP(ms)':>10}")
    print("-"*45)
    
    for n_waters in system_sizes:
        # 创建简单的测试系统
        state = pygcmc.MCState()
        atoms = []
        residues = []
        
        # 简单的立方体排列
        atoms_per_side = int(np.ceil(n_waters ** (1/3)))
        spacing = 0.3  # nm
        
        water_idx = 0
        for i in range(atoms_per_side):
            for j in range(atoms_per_side):
                for k in range(atoms_per_side):
                    if water_idx >= n_waters:
                        break
                    
                    x = i * spacing + 0.1
                    y = j * spacing + 0.1
                    z = k * spacing + 0.1
                    
                    # 添加水分子的5个原子
                    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
                    for atom_idx in range(5):
                        atom = pygcmc.MCAtom()
                        atom.x = x
                        atom.y = y
                        atom.z = z
                        atom.charge = charges[atom_idx]
                        atom.type = atom_idx if atom_idx < 4 else 3
                        atoms.append(atom)
                    
                    res = pygcmc.MCResidue()
                    res.atomStart = water_idx * 5
                    res.atomCount = 5
                    res.active = True
                    res.type = 0
                    residues.append(res)
                    
                    water_idx += 1
                if water_idx >= n_waters:
                    break
            if water_idx >= n_waters:
                break
        
        box_size = atoms_per_side * spacing + 0.2
        
        state.atoms = atoms
        state.residues = residues
        state.activeAtomCount = n_waters * 5
        state.activeResidueCount = n_waters
        
        state.info.box = np.array([box_size, box_size, box_size])
        state.info.cutoff = min(0.9, box_size / 2 - 0.01)
        
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 测试三种代表性算法
        times = {}
        
        for algo, name in [(pygcmc.DrudeAlgorithm.SCF, 'SCF'),
                          (pygcmc.DrudeAlgorithm.OPT3, 'OPT3'),
                          (pygcmc.DrudeAlgorithm.FBP, 'FBP')]:
            
            force = pygcmc.DrudeForce()
            
            for i in range(n_waters):
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
            
            force.setAlgorithm(algo)
            
            if algo == pygcmc.DrudeAlgorithm.SCF:
                params = pygcmc.DrudeSCFParams()
                params.tolerance = 100.0
                params.maxIterations = 20
                params.dampingFactor = 0.5
                params.maxDrudeDistance = 0.02
                force.setSCFParameters(params)
            
            state_test = state.copy()
            
            try:
                start = time.time()
                if algo == pygcmc.DrudeAlgorithm.SCF:
                    force.calculateEnergySCF(state_test)
                else:
                    force.calculateEnergyOPT(state_test)
                times[name] = (time.time() - start) * 1000
            except:
                times[name] = -1
        
        print(f"{n_waters:>10} {times.get('SCF', -1):>10.2f} {times.get('OPT3', -1):>10.2f} {times.get('FBP', -1):>10.2f}")

def main():
    """
    主函数
    """
    quick_benchmark()
    test_scaling()

if __name__ == "__main__":
    main()