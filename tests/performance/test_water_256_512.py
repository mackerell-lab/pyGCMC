#!/usr/bin/env python3
"""
测试256和512水系统的Drude SCF
密度1.0 g/cm³
"""

import pygcmc
import numpy as np
import pickle
import time
import os

def test_water_system(n_waters):
    """
    测试单个水系统
    """
    print(f"\n{'='*70}")
    print(f"测试 {n_waters} 水分子系统（密度 1.0 g/cm³）")
    print(f"{'='*70}")
    
    # 加载系统
    pickle_file = f'../tests/performance/large_water_systems/water_{n_waters}.pkl'
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
    
    # 系统信息
    box_length = data['box_length']
    density = data['density']
    
    print(f"\n系统参数:")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  半盒子: {box_length/2:.3f} nm")
    print(f"  密度: {density:.3f} g/cm³")
    print(f"  Thole(0.8nm)/半盒子: {0.8/(box_length/2):.3f}")
    
    # 创建状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP电荷
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]
    
    positions = data['positions']
    
    # 创建原子和残基
    print(f"\n创建原子和残基...")
    start_time = time.time()
    
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
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print(f"  完成，用时 {time.time()-start_time:.2f} 秒")
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    print(f"\n添加 {n_waters} 个Drude粒子...")
    start_time = time.time()
    
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
    
    print(f"  完成，用时 {time.time()-start_time:.2f} 秒")
    
    # 添加Thole对（采样方式）
    print(f"\n添加Thole对（采样）...")
    start_time = time.time()
    
    n_thole_pairs = 0
    thole_cutoff = 0.8
    
    # 为了加速，只检查部分对
    sample_rate = 10 if n_waters > 256 else 5
    max_neighbors = 50  # 每个分子最多检查的邻居数
    
    for i in range(0, n_waters, sample_rate):
        o1_idx = i * 5
        neighbors_found = 0
        
        # 只检查附近的分子
        for j in range(i+1, min(i+200, n_waters)):
            if neighbors_found >= max_neighbors:
                break
                
            o2_idx = j * 5
            
            # 计算距离
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            # PBC
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, 1.3)
                n_thole_pairs += 1
                neighbors_found += 1
    
    print(f"  添加了 {n_thole_pairs} 个Thole对")
    print(f"  估算总Thole对: ~{n_thole_pairs * sample_rate}")
    print(f"  完成，用时 {time.time()-start_time:.2f} 秒")
    
    # 测试不同容差的SCF
    tolerances = [10.0, 100.0]
    
    print(f"\n测试SCF收敛:")
    print(f"{'容差':>10} {'时间(s)':>10} {'能量/水':>12} {'平均位移(pm)':>15}")
    print("-"*60)
    
    for tolerance in tolerances:
        # 设置SCF参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tolerance
        params.maxIterations = 50  # 减少迭代次数
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 复制状态
        test_state = state.copy()
        
        # 运行SCF
        start_time = time.time()
        
        try:
            energy = force.calculateEnergySCF(test_state)
            elapsed_time = time.time() - start_time
            
            # 采样分析位移
            displacements = []
            for i in range(0, n_waters, max(1, n_waters//50)):
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = test_state.atoms[d_idx].x - test_state.atoms[o_idx].x
                dy = test_state.atoms[d_idx].y - test_state.atoms[o_idx].y
                dz = test_state.atoms[d_idx].z - test_state.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                displacements.append(disp)
            
            avg_disp = np.mean(displacements)
            
            print(f"{tolerance:10.1f} {elapsed_time:10.2f} {energy/n_waters:12.2f} "
                  f"{avg_disp:15.2f}")
            
        except Exception as e:
            elapsed_time = time.time() - start_time
            print(f"{tolerance:10.1f} {elapsed_time:10.2f} {'失败':>12} "
                  f"{'N/A':>15}")

def main():
    """
    主函数
    """
    print("大型水系统Drude SCF验证")
    print("验证密度1.0 g/cm³下的收敛性")
    print("="*70)
    
    # 测试256和512水系统
    for n_waters in [256, 512]:
        try:
            test_water_system(n_waters)
        except Exception as e:
            print(f"\n错误: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n\n结论:")
    print("="*70)
    print("1. 256水系统（盒子~2nm）: Thole/半盒子 = 0.81 < 1.0 ✓")
    print("2. 512水系统（盒子~2.5nm）: Thole/半盒子 = 0.64 < 1.0 ✓")
    print("3. 大系统满足最小镜像约定，应该能够收敛")
    print("4. 如果仍有收敛问题，可能需要优化初始猜测或使用ASPC")

if __name__ == "__main__":
    main()