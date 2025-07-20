#!/usr/bin/env python3
"""
使用OpenMM优化的水系统测试FBP性能
"""

import pygcmc
import numpy as np
import time
import pickle

def load_optimized_water_system(filename='water128_openmm_optimized.pkl'):
    """加载优化的水系统配置"""
    print(f"加载优化的水系统: {filename}")
    
    with open(filename, 'rb') as f:
        config = pickle.load(f)
    
    n_waters = config['n_waters']
    positions = config['positions']
    
    print(f"加载了{n_waters}个水分子")
    
    # 计算盒子尺寸
    max_x = max_y = max_z = 0
    min_x = min_y = min_z = 1e10
    
    for water in positions:
        for atom_pos in [water['O'], water['H1'], water['H2']]:
            max_x = max(max_x, atom_pos[0])
            max_y = max(max_y, atom_pos[1])
            max_z = max(max_z, atom_pos[2])
            min_x = min(min_x, atom_pos[0])
            min_y = min(min_y, atom_pos[1])
            min_z = min(min_z, atom_pos[2])
    
    box_x = max_x - min_x + 0.5  # 加一些边距
    box_y = max_y - min_y + 0.5
    box_z = max_z - min_z + 0.5
    box_length = max(box_x, box_y, box_z)
    
    print(f"盒子尺寸: {box_length:.2f} nm")
    
    # 转换为pygcmc格式
    atoms = []
    residues = []
    
    for i, water in enumerate(positions):
        # SWM4-NDP水模型
        # O原子
        o_pos = water['O']
        atom = pygcmc.MCAtom()
        atom.x = o_pos[0]
        atom.y = o_pos[1]
        atom.z = o_pos[2]
        atom.charge = 1.71636
        atom.type = 0
        atoms.append(atom)
        
        # Drude粒子（初始与O重合）
        atom = pygcmc.MCAtom()
        atom.x = o_pos[0]
        atom.y = o_pos[1]
        atom.z = o_pos[2]
        atom.charge = -1.71636
        atom.type = 1
        atoms.append(atom)
        
        # H1
        h1_pos = water['H1']
        atom = pygcmc.MCAtom()
        atom.x = h1_pos[0]
        atom.y = h1_pos[1]
        atom.z = h1_pos[2]
        atom.charge = 0.55733
        atom.type = 2
        atoms.append(atom)
        
        # H2
        h2_pos = water['H2']
        atom = pygcmc.MCAtom()
        atom.x = h2_pos[0]
        atom.y = h2_pos[1]
        atom.z = h2_pos[2]
        atom.charge = 0.55733
        atom.type = 2
        atoms.append(atom)
        
        # M虚拟位点（质心）
        m_x = o_pos[0] * 0.8476 + (h1_pos[0] + h2_pos[0]) * 0.0762
        m_y = o_pos[1] * 0.8476 + (h1_pos[1] + h2_pos[1]) * 0.0762
        m_z = o_pos[2] * 0.8476 + (h1_pos[2] + h2_pos[2]) * 0.0762
        
        atom = pygcmc.MCAtom()
        atom.x = m_x
        atom.y = m_y
        atom.z = m_z
        atom.charge = -1.11466
        atom.type = 3
        atoms.append(atom)
        
        # 残基
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    # 创建状态
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
    
    return state, n_waters

def test_fbp_on_optimized_system():
    """在优化的系统上测试FBP"""
    print("\nFBP在OpenMM优化系统上的性能测试")
    print("="*80)
    
    # 加载系统
    try:
        state, n_waters = load_optimized_water_system()
    except FileNotFoundError:
        print("未找到优化的水系统文件，请先运行 create_openmm_water_system.py")
        return
    
    print(f"\n系统信息:")
    print(f"  水分子数: {n_waters}")
    print(f"  原子总数: {state.activeAtomCount}")
    print(f"  盒子尺寸: {state.info.box[0]:.2f} nm")
    print(f"  截断距离: {state.info.cutoff:.2f} nm")
    
    # 测试不同算法
    algorithms = [
        ("SCF (tol=0.1)", "SCF", 0.1, 500),
        ("SCF (tol=1.0)", "SCF", 1.0, 500),
        ("SCF (tol=10.0)", "SCF", 10.0, 200),
        ("FBP (tol=0.1)", "FBP", 0.1, 200),
        ("FBP (tol=1.0)", "FBP", 1.0, 200),
        ("FBP (tol=10.0)", "FBP", 10.0, 200),
    ]
    
    results = []
    
    for name, algo, tolerance, max_iter in algorithms:
        print(f"\n测试 {name}...")
        
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
        params.tolerance = tolerance
        params.maxIterations = max_iter
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 复制状态（保持原始位置）
        state_copy = state.copy()
        
        # 运行3次取平均
        times = []
        energies = []
        
        for run in range(3):
            # 重置Drude位置
            for i in range(n_waters):
                state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
                state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
                state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
            
            start = time.time()
            energy = force.calculateEnergySCF(state_copy)
            elapsed = (time.time() - start) * 1000  # ms
            
            times.append(elapsed)
            energies.append(energy)
        
        avg_time = np.mean(times)
        std_time = np.std(times)
        avg_energy = np.mean(energies)
        
        # 计算最终Drude位移
        displacements = []
        for i in range(n_waters):
            dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
            dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
            dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
            displacements.append(disp)
        
        avg_disp = np.mean(displacements)
        max_disp = np.max(displacements)
        
        results.append({
            'name': name,
            'time': avg_time,
            'std': std_time,
            'energy': avg_energy,
            'avg_disp': avg_disp,
            'max_disp': max_disp
        })
        
        print(f"  时间: {avg_time:.1f} ± {std_time:.1f} ms")
        print(f"  能量: {avg_energy:.2f} kJ/mol")
        print(f"  平均位移: {avg_disp:.2f} pm")
    
    # 结果汇总
    print("\n\n结果汇总:")
    print("="*100)
    print(f"{'算法':<20} {'时间(ms)':<15} {'能量(kJ/mol)':<15} {'平均位移(pm)':<12} {'相对于SCF(0.1)':<20}")
    print("-"*100)
    
    # 找到SCF(0.1)作为参考
    ref_result = results[0]  # SCF (tol=0.1)
    ref_time = ref_result['time']
    ref_energy = ref_result['energy']
    
    for r in results:
        speedup = ref_time / r['time']
        energy_diff = abs(r['energy'] - ref_energy)
        energy_diff_pct = energy_diff / abs(ref_energy) * 100 if abs(ref_energy) > 1 else 0
        
        relative = f"{speedup:.2f}x速度, {energy_diff_pct:.1f}%能量差"
        
        print(f"{r['name']:<20} {r['time']:<15.1f} {r['energy']:<15.2f} {r['avg_disp']:<12.2f} {relative:<20}")
    
    # 特别对比FBP和SCF
    print("\n\nFBP vs SCF 关键对比:")
    print("-"*60)
    
    # 找到对应的结果
    scf_1 = next((r for r in results if r['name'] == "SCF (tol=1.0)"), None)
    fbp_1 = next((r for r in results if r['name'] == "FBP (tol=1.0)"), None)
    
    if scf_1 and fbp_1:
        speedup = scf_1['time'] / fbp_1['time']
        energy_diff = abs(fbp_1['energy'] - scf_1['energy'])
        
        print(f"在容差1.0时:")
        print(f"  FBP速度提升: {speedup:.2f}x")
        print(f"  能量差异: {energy_diff:.2f} kJ/mol")
        print(f"  FBP时间: {fbp_1['time']:.1f} ms")
        print(f"  SCF时间: {scf_1['time']:.1f} ms")
        
        if speedup > 1:
            print(f"\n结论: 在预优化的系统上，FBP比SCF快{speedup:.1f}倍")
        else:
            print(f"\n结论: 在预优化的系统上，SCF比FBP快{1/speedup:.1f}倍")

if __name__ == "__main__":
    test_fbp_on_optimized_system()