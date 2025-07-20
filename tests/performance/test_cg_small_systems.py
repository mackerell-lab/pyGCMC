#!/usr/bin/env python3
"""
CG vs SCF对比测试 - 仅测试小系统
"""

import pygcmc
import numpy as np
import time
import pickle
import os

def load_water_system(n_waters):
    """
    加载水体系
    """
    pickle_file = f'../tests/performance/water_systems/water_{n_waters}.pkl'
    if os.path.exists(pickle_file):
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
        return convert_to_pygcmc(data)
    else:
        raise FileNotFoundError(f"找不到: {pickle_file}")

def convert_to_pygcmc(data):
    """
    将保存的数据转换为pygcmc状态
    """
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters):
        for j in range(5):
            atom = pygcmc.MCAtom()
            idx = i * 5 + j
            atom.x = positions[idx][0]
            atom.y = positions[idx][1]
            atom.z = positions[idx][2]
            atom.charge = charges[j]
            atom.type = atom_types[j]
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
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length / 2 - 0.01)
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_single_system(n_waters):
    """
    测试单个系统
    """
    print(f"\n测试 {n_waters} 水分子系统")
    print("-"*50)
    
    # 加载系统
    state = load_water_system(n_waters)
    
    # 系统信息
    actual_density = n_waters * 18.015 / (state.info.box[0]**3 * 0.6022)
    print(f"盒子: {state.info.box[0]:.3f} nm, 密度: {actual_density:.3f} g/cm³")
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    polarizability = 0.0009782237  # nm³
    
    # 添加Drude粒子
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 添加Thole对
    n_pairs = 0
    cutoff_thole = 0.8  # nm
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            dx = state.atoms[5*i].x - state.atoms[5*j].x
            dy = state.atoms[5*i].y - state.atoms[5*j].y
            dz = state.atoms[5*i].z - state.atoms[5*j].z
            
            # PBC
            box = state.info.box[0]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < cutoff_thole:
                force.addScreenedPair(i, j, 1.3)
                n_pairs += 1
    
    print(f"Thole对: {n_pairs}")
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # SCF测试
    print("\nSCF:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    state_scf = state.copy()
    
    start = time.time()
    energy_scf = force.calculateEnergySCF(state_scf)
    time_scf = (time.time() - start) * 1000
    
    # 分析位移
    disps_scf = []
    for i in range(n_waters):
        dx = state_scf.atoms[5*i+1].x - state_scf.atoms[5*i].x
        dy = state_scf.atoms[5*i+1].y - state_scf.atoms[5*i].y
        dz = state_scf.atoms[5*i+1].z - state_scf.atoms[5*i].z
        disps_scf.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
    
    print(f"  时间: {time_scf:.1f} ms")
    print(f"  能量: {energy_scf:.2f} kJ/mol ({energy_scf/n_waters:.2f} per water)")
    print(f"  平均位移: {np.mean(disps_scf):.2f} pm")
    
    # CG测试
    print("\nCG:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
    state_cg = state.copy()
    
    start = time.time()
    energy_cg = force.calculateEnergySCF(state_cg)
    time_cg = (time.time() - start) * 1000
    
    # 分析位移
    disps_cg = []
    for i in range(n_waters):
        dx = state_cg.atoms[5*i+1].x - state_cg.atoms[5*i].x
        dy = state_cg.atoms[5*i+1].y - state_cg.atoms[5*i].y
        dz = state_cg.atoms[5*i+1].z - state_cg.atoms[5*i].z
        disps_cg.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
    
    print(f"  时间: {time_cg:.1f} ms")
    print(f"  能量: {energy_cg:.2f} kJ/mol ({energy_cg/n_waters:.2f} per water)")
    print(f"  平均位移: {np.mean(disps_cg):.2f} pm")
    
    # 对比
    print("\n对比:")
    speedup = time_scf / time_cg if time_cg > 0 else 0
    energy_diff = abs(energy_cg - energy_scf)
    energy_diff_pct = energy_diff / abs(energy_scf) * 100 if energy_scf != 0 else 0
    
    print(f"  加速: {speedup:.2f}x")
    print(f"  能量差: {energy_diff:.2f} kJ/mol ({energy_diff_pct:.1f}%)")
    print(f"  位移差: {abs(np.mean(disps_cg) - np.mean(disps_scf)):.2f} pm")

def main():
    """
    主函数
    """
    print("CG vs SCF 小系统测试")
    print("="*50)
    
    # 只测试小系统
    system_sizes = [2, 4, 8, 16, 32]
    
    for n_waters in system_sizes:
        try:
            test_single_system(n_waters)
        except Exception as e:
            print(f"\n错误: {n_waters} 水失败: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    main()