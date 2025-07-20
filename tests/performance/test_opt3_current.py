#!/usr/bin/env python3
"""
测试OPT3算法的当前性能
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
    # 尝试不同的目录
    for directory in ['../tests/performance/water_systems_final',
                     '../tests/performance/water_systems']:
        pickle_file = os.path.join(directory, f'water_{n_waters}.pkl')
        if os.path.exists(pickle_file):
            with open(pickle_file, 'rb') as f:
                data = pickle.load(f)
            return convert_to_pygcmc(data)
    raise FileNotFoundError(f"找不到水系统文件")

def convert_to_pygcmc(data):
    """
    转换为pygcmc格式
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

def test_opt3_performance(n_waters):
    """
    测试OPT3性能
    """
    print(f"\n测试 {n_waters} 水分子系统")
    print("-"*60)
    
    # 加载系统
    state = load_water_system(n_waters)
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    charge = -1.71636
    polarizability = 0.0009782237  # nm³
    
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
    
    print(f"系统信息:")
    print(f"  盒子: {state.info.box[0]:.3f} nm")
    print(f"  密度: {n_waters * 18.015 / (state.info.box[0]**3 * 0.6022):.3f} g/cm³")
    print(f"  Thole对: {n_pairs}")
    
    # 获取当前OPT3系数
    coeffs = force.getOPT3Coefficients()
    print(f"\n当前OPT3系数:")
    print(f"  c0 = {coeffs.c0}")
    print(f"  c1 = {coeffs.c1}")
    print(f"  c2 = {coeffs.c2}")
    print(f"  c3 = {coeffs.c3}")
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 测试SCF
    print("\nSCF:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    state_scf = state.copy()
    
    start = time.time()
    energy_scf = force.calculateEnergySCF(state_scf)
    time_scf = (time.time() - start) * 1000
    
    print(f"  时间: {time_scf:.1f} ms")
    print(f"  能量: {energy_scf:.2f} kJ/mol")
    
    # 测试OPT3 (默认系数)
    print("\nOPT3 (默认系数):")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    state_opt3_default = state.copy()
    
    start = time.time()
    energy_opt3_default = force.calculateEnergySCF(state_opt3_default)
    time_opt3_default = (time.time() - start) * 1000
    
    print(f"  时间: {time_opt3_default:.1f} ms")
    print(f"  能量: {energy_opt3_default:.2f} kJ/mol")
    print(f"  能量差: {abs(energy_opt3_default - energy_scf):.2f} kJ/mol ({abs(energy_opt3_default - energy_scf)/abs(energy_scf)*100:.1f}%)")
    print(f"  速度: {time_scf/time_opt3_default:.2f}x")
    
    # 测试优化的OPT3系数
    print("\nOPT3 (优化系数 c0=2.0, c1=1.0, c2=0.0, c3=0.0):")
    force.setOPT3Coefficients(2.0, 1.0, 0.0, 0.0)
    state_opt3_optimized = state.copy()
    
    start = time.time()
    energy_opt3_optimized = force.calculateEnergySCF(state_opt3_optimized)
    time_opt3_optimized = (time.time() - start) * 1000
    
    print(f"  时间: {time_opt3_optimized:.1f} ms")
    print(f"  能量: {energy_opt3_optimized:.2f} kJ/mol")
    print(f"  能量差: {abs(energy_opt3_optimized - energy_scf):.2f} kJ/mol ({abs(energy_opt3_optimized - energy_scf)/abs(energy_scf)*100:.1f}%)")
    print(f"  速度: {time_scf/time_opt3_optimized:.2f}x")

def main():
    """
    主函数
    """
    print("OPT3算法性能测试")
    print("="*60)
    
    # 测试不同大小的系统
    system_sizes = [2, 4, 8, 16, 32]
    
    for n_waters in system_sizes:
        try:
            test_opt3_performance(n_waters)
        except Exception as e:
            print(f"\n错误: {n_waters} 水分子失败: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    main()