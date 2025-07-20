#!/usr/bin/env python3
"""
测试FBP（Force Balance Propagation）方法
"""

import numpy as np
import time
import pickle
import os
import pygcmc

def test_fbp_method():
    """
    测试FBP方法的正确调用方式
    """
    print("测试FBP方法")
    print("="*70)
    
    # 创建简单的2水系统
    state = pygcmc.MCState()
    
    box_size = 1.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
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
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
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
    
    force.addScreenedPair(0, 1, 1.3)
    
    # 测试不同的方法名
    print("\n1. 查找FBP相关方法:")
    methods = [m for m in dir(force) if not m.startswith('_')]
    fbp_methods = [m for m in methods if 'fbp' in m.lower() or 'FBP' in m]
    print(f"   FBP相关方法: {fbp_methods}")
    
    energy_methods = [m for m in methods if 'calculateEnergy' in m]
    print(f"   所有能量计算方法: {energy_methods}")
    
    # 测试FBP算法
    print("\n2. 测试FBP算法:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 尝试不同的方法
    test_methods = ['calculateEnergyFBP', 'calculateEnergyOPT', 'calculateEnergyOPT3', 'calculateEnergy']
    
    for method_name in test_methods:
        if hasattr(force, method_name):
            print(f"\n   尝试 {method_name}...")
            try:
                state_test = state.copy()
                method = getattr(force, method_name)
                energy = method(state_test)
                
                # 计算位移
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
                
                print(f"   ✓ 成功! 能量={energy:.4f} kJ/mol, 平均位移={avg_disp:.2f} pm")
                break
            except Exception as e:
                print(f"   ✗ 失败: {str(e)[:50]}")
        else:
            print(f"   - {method_name} 不存在")

def benchmark_all_methods():
    """
    基准测试所有可用方法
    """
    print("\n\n完整性能测试（SCF vs OPT3 vs FBP）")
    print("="*70)
    
    # 加载测试系统
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    if not os.path.exists(filename):
        print("文件不存在")
        return
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    # 测试不同大小
    test_sizes = [10, 50, 100, 256]
    
    print(f"{'系统大小':>10} {'算法':>10} {'时间(ms)':>12} {'能量/水':>15} {'位移(pm)':>12} {'速度比':>10}")
    print("-"*80)
    
    for n_waters in test_sizes:
        if n_waters > 256:
            continue
            
        # 创建状态
        state = create_state_from_data(data, n_waters)
        
        # 测试三种算法
        results = {}
        
        algorithms = [
            (pygcmc.DrudeAlgorithm.SCF, "SCF", "calculateEnergySCF"),
            (pygcmc.DrudeAlgorithm.OPT3, "OPT3", "calculateEnergyOPT3"),
            (pygcmc.DrudeAlgorithm.FBP, "FBP", None)  # 待确定
        ]
        
        for algo, algo_name, method_name in algorithms:
            if algo_name == "FBP" and method_name is None:
                # 尝试找到正确的FBP方法
                continue
                
            force = create_drude_force(state, n_waters, data['box_length'])
            force.setAlgorithm(algo)
            
            if algo == pygcmc.DrudeAlgorithm.SCF:
                params = pygcmc.DrudeSCFParams()
                params.tolerance = 100.0
                params.maxIterations = 30
                params.dampingFactor = 0.5
                params.maxDrudeDistance = 0.02
                force.setSCFParameters(params)
            
            # 准备状态
            state_test = state.copy()
            reset_drude_positions(state_test, n_waters)
            
            try:
                # 计时
                start = time.time()
                method = getattr(force, method_name)
                energy = method(state_test)
                elapsed = (time.time() - start) * 1000
                
                # 计算位移
                avg_disp = calculate_avg_displacement(state_test, n_waters)
                
                results[algo_name] = {
                    'time': elapsed,
                    'energy': energy / n_waters,
                    'displacement': avg_disp
                }
                
            except Exception as e:
                results[algo_name] = {'error': str(e)}
        
        # 打印结果
        if 'SCF' in results and 'time' in results['SCF']:
            scf_time = results['SCF']['time']
        else:
            scf_time = 1.0
            
        for algo_name in ['SCF', 'OPT3', 'FBP']:
            if algo_name in results:
                r = results[algo_name]
                if 'error' in r:
                    print(f"{n_waters:>10} {algo_name:>10} {'失败':>12} {'-':>15} {'-':>12} {'-':>10}")
                else:
                    speed_ratio = r['time'] / scf_time
                    print(f"{n_waters:>10} {algo_name:>10} {r['time']:>12.2f} {r['energy']:>15.2f} {r['displacement']:>12.2f} {speed_ratio:>10.2f}x")

def create_state_from_data(data, n_waters):
    """创建状态"""
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    positions = data['positions']
    charges = data['charges']
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
    
    state.info.box = np.array([data['box_length']] * 3)
    state.info.cutoff = min(0.9, data['box_length'] / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def create_drude_force(state, n_waters, box_length):
    """创建DrudeForce"""
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
    
    # 添加Thole对
    for i in range(n_waters):
        for j in range(i+1, min(i+10, n_waters)):
            o1_idx = i * 5
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < 0.8:
                force.addScreenedPair(i, j, 1.3)
    
    return force

def reset_drude_positions(state, n_waters):
    """重置Drude位置"""
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        state.atoms[d_idx].x = state.atoms[o_idx].x
        state.atoms[d_idx].y = state.atoms[o_idx].y
        state.atoms[d_idx].z = state.atoms[o_idx].z

def calculate_avg_displacement(state, n_waters):
    """计算平均位移"""
    displacements = []
    for i in range(min(20, n_waters)):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = state.atoms[d_idx].x - state.atoms[o_idx].x
        dy = state.atoms[d_idx].y - state.atoms[o_idx].y
        dz = state.atoms[d_idx].z - state.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        displacements.append(disp)
    
    return np.mean(displacements)

def main():
    """
    主函数
    """
    test_fbp_method()
    benchmark_all_methods()
    
    print("\n\n结论：")
    print("="*70)
    print("需要确定FBP方法的正确调用方式")

if __name__ == "__main__":
    main()