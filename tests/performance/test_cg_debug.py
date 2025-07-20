#!/usr/bin/env python3
"""
调试CG实现
"""

import pygcmc
import numpy as np

def create_single_water():
    """创建单个水分子"""
    atoms = []
    
    # SWM4-NDP水模型
    positions = [
        (0.5, 0.5, 0.5, 1.71636, 0),   # O
        (0.5, 0.5, 0.5, -1.71636, 1),  # D (初始与O重合)
        (0.59572, 0.5, 0.5, 0.55733, 2),  # H1
        (0.45214, 0.58288, 0.5, 0.55733, 2),  # H2
        (0.5, 0.475966, 0.5, -1.11466, 3)  # M
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
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = [res]
    state.activeAtomCount = 5
    state.activeResidueCount = 1
    
    state.info.box = np.array([2.0, 2.0, 2.0])
    state.info.cutoff = 0.9
    
    # SWM4-NDP力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_single_water():
    """测试单个水分子"""
    state = create_single_water()
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    print(f"Drude参数:")
    print(f"  charge = {charge}")
    print(f"  k_spring = {k_spring}")
    print(f"  polarizability = {polarizability}")
    
    # 添加Drude粒子
    force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 100.0  # 使用很大容差
    params.maxIterations = 10  # 很少迭代
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1  # 增大硬墙距离
    force.setSCFParameters(params)
    
    # 测试SCF
    print("\n测试SCF:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    state_scf = state.copy()
    energy_scf = force.calculateEnergySCF(state_scf)
    
    dx_scf = state_scf.atoms[1].x - state_scf.atoms[0].x
    dy_scf = state_scf.atoms[1].y - state_scf.atoms[0].y
    dz_scf = state_scf.atoms[1].z - state_scf.atoms[0].z
    disp_scf = np.sqrt(dx_scf**2 + dy_scf**2 + dz_scf**2) * 1000
    
    print(f"  能量: {energy_scf:.2f} kJ/mol")
    print(f"  位移: {disp_scf:.2f} pm")
    print(f"  dx,dy,dz: {dx_scf*1000:.3f}, {dy_scf*1000:.3f}, {dz_scf*1000:.3f} pm")
    
    # 测试CG
    print("\n测试CG:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
    state_cg = state.copy()
    energy_cg = force.calculateEnergySCF(state_cg)
    
    dx_cg = state_cg.atoms[1].x - state_cg.atoms[0].x
    dy_cg = state_cg.atoms[1].y - state_cg.atoms[0].y
    dz_cg = state_cg.atoms[1].z - state_cg.atoms[0].z
    disp_cg = np.sqrt(dx_cg**2 + dy_cg**2 + dz_cg**2) * 1000
    
    print(f"  能量: {energy_cg:.2f} kJ/mol")
    print(f"  位移: {disp_cg:.2f} pm")
    print(f"  dx,dy,dz: {dx_cg*1000:.3f}, {dy_cg*1000:.3f}, {dz_cg*1000:.3f} pm")
    
    print(f"\n差异:")
    print(f"  能量差: {abs(energy_cg - energy_scf):.2f} kJ/mol")
    print(f"  位移差: {abs(disp_cg - disp_scf):.2f} pm")

def test_two_waters():
    """测试两个水分子相互作用"""
    atoms = []
    
    # 第一个水分子
    positions1 = [
        (0.5, 0.5, 0.5, 1.71636, 0),   # O
        (0.5, 0.5, 0.5, -1.71636, 1),  # D
        (0.59572, 0.5, 0.5, 0.55733, 2),  # H1
        (0.45214, 0.58288, 0.5, 0.55733, 2),  # H2
        (0.5, 0.475966, 0.5, -1.11466, 3)  # M
    ]
    
    # 第二个水分子
    positions2 = [
        (0.8, 0.5, 0.5, 1.71636, 0),   # O
        (0.8, 0.5, 0.5, -1.71636, 1),  # D
        (0.89572, 0.5, 0.5, 0.55733, 2),  # H1
        (0.75214, 0.58288, 0.5, 0.55733, 2),  # H2
        (0.8, 0.475966, 0.5, -1.11466, 3)  # M
    ]
    
    for positions in [positions1, positions2]:
        for px, py, pz, charge, typ in positions:
            atom = pygcmc.MCAtom()
            atom.x = px
            atom.y = py
            atom.z = pz
            atom.charge = charge
            atom.type = typ
            atoms.append(atom)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeAtomCount = 10
    state.activeResidueCount = 2
    
    state.info.box = np.array([2.0, 2.0, 2.0])
    state.info.cutoff = 1.5
    
    # SWM4-NDP力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    # 添加两个Drude粒子
    for i in range(2):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 添加Thole屏蔽
    force.addScreenedPair(0, 1, 1.3)
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 100.0
    params.maxIterations = 20
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1
    force.setSCFParameters(params)
    
    print("\n\n测试两个水分子:")
    
    # 测试SCF
    print("\nSCF:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    state_scf = state.copy()
    energy_scf = force.calculateEnergySCF(state_scf)
    print(f"  能量: {energy_scf:.2f} kJ/mol")
    
    # 测试CG
    print("\nCG:")
    force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
    state_cg = state.copy()
    energy_cg = force.calculateEnergySCF(state_cg)
    print(f"  能量: {energy_cg:.2f} kJ/mol")
    
    print(f"\n能量差: {abs(energy_cg - energy_scf):.2f} kJ/mol")

if __name__ == "__main__":
    test_single_water()
    test_two_waters()