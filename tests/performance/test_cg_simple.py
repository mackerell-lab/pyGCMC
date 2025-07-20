#!/usr/bin/env python3
"""
简单测试共轭梯度法是否正常工作
"""

import pygcmc
import numpy as np
import time

def create_simple_water_system(n_waters=4):
    """创建简单水分子系统"""
    spacing = 0.5  # nm
    
    atoms = []
    residues = []
    
    water_count = 0
    for i in range(n_waters):
        x = (i + 0.5) * spacing
        y = 0.5
        z = 0.5
        
        # SWM4-NDP水模型
        positions = [
            (x, y, z, 1.71636, 0),   # O
            (x, y, z, -1.71636, 1),  # D (初始与O重合)
            (x + 0.09572, y, z, 0.55733, 2),  # H1
            (x - 0.04786, y + 0.08288, z, 0.55733, 2),  # H2
            (x, y - 0.024034, z, -1.11466, 3)  # M
        ]
        
        for px, py, pz, charge, typ in positions:
            atom = pygcmc.MCAtom()
            atom.x = px
            atom.y = py
            atom.z = pz
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
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    box_length = n_waters * spacing + 1.0
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length/2 - 0.01)
    
    # SWM4-NDP力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_algorithm(state, n_waters, algorithm):
    """测试算法"""
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
    
    # 添加Thole屏蔽（只添加少量）
    if n_waters > 1:
        force.addScreenedPair(0, 1, 1.3)
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # 使用较大容差
    params.maxIterations = 50  # 减少迭代次数
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 设置算法
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm == "CG":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    
    print(f"\nTesting {algorithm} with {n_waters} waters...")
    
    # 重置Drude位置
    for i in range(n_waters):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    start = time.time()
    energy = force.calculateEnergySCF(state)
    elapsed = (time.time() - start) * 1000  # ms
    
    # 计算最终位移
    displacements = []
    for i in range(n_waters):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        dy = state.atoms[5*i+1].y - state.atoms[5*i].y
        dz = state.atoms[5*i+1].z - state.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        displacements.append(disp)
    
    avg_disp = np.mean(displacements)
    
    print(f"  时间: {elapsed:.1f} ms")
    print(f"  能量: {energy:.2f} kJ/mol")
    print(f"  平均位移: {avg_disp:.2f} pm")
    
    return elapsed, energy, avg_disp

def main():
    """主测试函数"""
    print("共轭梯度法(CG)简单测试")
    print("="*50)
    
    # 测试不同大小的系统
    for n_waters in [2, 4, 8]:
        state = create_simple_water_system(n_waters)
        
        # 测试两种算法
        scf_time, scf_energy, scf_disp = test_algorithm(state, n_waters, "SCF")
        cg_time, cg_energy, cg_disp = test_algorithm(state, n_waters, "CG")
        
        print(f"\n{n_waters}水分子系统总结:")
        print(f"  SCF: {scf_time:.1f} ms, 能量 = {scf_energy:.2f} kJ/mol")
        print(f"  CG:  {cg_time:.1f} ms, 能量 = {cg_energy:.2f} kJ/mol")
        print(f"  加速比: {scf_time/cg_time:.2f}x")
        print(f"  能量差: {abs(cg_energy - scf_energy):.2f} kJ/mol")

if __name__ == "__main__":
    main()