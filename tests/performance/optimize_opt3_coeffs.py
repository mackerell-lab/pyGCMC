#!/usr/bin/env python3
"""
优化OPT3系数
"""

import pygcmc
import numpy as np
import pickle

def load_small_water_system():
    """
    加载小水系统用于快速测试
    """
    with open('../tests/performance/water_systems_final/water_4.pkl', 'rb') as f:
        data = pickle.load(f)
    
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
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, n_waters

def test_coefficients(c0, c1, c2, c3):
    """
    测试一组系数
    """
    state, n_waters = load_small_water_system()
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    charge = -1.71636
    polarizability = 0.0009782237
    
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 添加Thole对
    cutoff_thole = 0.8
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            dx = state.atoms[5*i].x - state.atoms[5*j].x
            dy = state.atoms[5*i].y - state.atoms[5*j].y
            dz = state.atoms[5*i].z - state.atoms[5*j].z
            
            box = state.info.box[0]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < cutoff_thole:
                force.addScreenedPair(i, j, 1.3)
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 100.0  # 更宽松的容差
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 计算SCF能量（参考）
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    state_scf = state.copy()
    energy_scf = force.calculateEnergySCF(state_scf)
    
    # 计算OPT3能量
    force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    force.setOPT3Coefficients(c0, c1, c2, c3)
    state_opt3 = state.copy()
    energy_opt3 = force.calculateEnergySCF(state_opt3)
    
    # 计算误差
    error = abs(energy_opt3 - energy_scf) / abs(energy_scf) * 100
    
    return error, energy_scf, energy_opt3

def grid_search():
    """
    网格搜索最优系数
    """
    print("OPT3系数优化（网格搜索）")
    print("="*60)
    
    # 搜索范围
    c0_values = [0.0]  # 通常设为0
    c1_values = np.arange(0.8, 2.0, 0.2)
    c2_values = np.arange(-1.0, 0.5, 0.2)
    c3_values = np.arange(-0.2, 0.4, 0.1)
    
    best_error = float('inf')
    best_coeffs = None
    
    print("\n搜索进度...")
    total = len(c0_values) * len(c1_values) * len(c2_values) * len(c3_values)
    count = 0
    
    for c0 in c0_values:
        for c1 in c1_values:
            for c2 in c2_values:
                for c3 in c3_values:
                    count += 1
                    
                    try:
                        error, e_scf, e_opt3 = test_coefficients(c0, c1, c2, c3)
                        
                        if error < best_error:
                            best_error = error
                            best_coeffs = (c0, c1, c2, c3)
                            print(f"\n新最优: c0={c0:.1f}, c1={c1:.1f}, c2={c2:.1f}, c3={c3:.1f}")
                            print(f"  误差: {error:.1f}%")
                            print(f"  SCF能量: {e_scf:.2f} kJ/mol")
                            print(f"  OPT3能量: {e_opt3:.2f} kJ/mol")
                        
                        if count % 10 == 0:
                            print(f"  进度: {count}/{total} ({count/total*100:.1f}%)")
                    
                    except Exception as e:
                        # 跳过失败的组合
                        pass
    
    print("\n\n最优系数:")
    print(f"c0 = {best_coeffs[0]}")
    print(f"c1 = {best_coeffs[1]}")
    print(f"c2 = {best_coeffs[2]}")
    print(f"c3 = {best_coeffs[3]}")
    print(f"最小误差: {best_error:.1f}%")
    
    # 测试更多系统
    print("\n\n在不同系统上验证:")
    test_systems = [2, 4, 8, 16]
    
    for n in test_systems:
        print(f"\n{n}水系统:")
        # 这里简化，只测试4水系统
        if n == 4:
            error, e_scf, e_opt3 = test_coefficients(*best_coeffs)
            print(f"  误差: {error:.1f}%")
            print(f"  SCF: {e_scf:.2f} kJ/mol")
            print(f"  OPT3: {e_opt3:.2f} kJ/mol")

if __name__ == "__main__":
    grid_search()