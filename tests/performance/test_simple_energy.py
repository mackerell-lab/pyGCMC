#!/usr/bin/env python3
"""
测试简单系统的能量计算
"""

import pygcmc
import numpy as np

def test_simple_energy():
    """测试最简单情况的能量"""
    print("简单能量测试")
    print("="*60)
    
    # 只测试两个水分子，不同间距
    distances = [0.3, 0.4, 0.5, 0.7, 1.0, 1.5]
    
    print(f"\n{'间距(nm)':<10} {'SCF能量':<15} {'FBP能量':<15} {'差异':<10}")
    print("-"*50)
    
    for dist in distances:
        # 创建系统
        atoms = []
        residues = []
        
        # 两个水分子
        for i in range(2):
            x_base = i * dist
            positions = [
                (x_base, 0.0, 0.0, 1.71636, 0),   # O
                (x_base, 0.0, 0.0, -1.71636, 1),  # D
                (x_base + 0.09572, 0.0, 0.0, 0.55733, 2),  # H1
                (x_base - 0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
                (x_base, -0.024034, 0.0, -1.11466, 3)  # M
            ]
            
            for x, y, z, charge, typ in positions:
                atom = pygcmc.MCAtom()
                atom.x = x
                atom.y = y
                atom.z = z
                atom.charge = charge
                atom.type = typ
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
        state.activeAtomCount = 10
        state.activeResidueCount = 2
        
        state.info.box = np.array([10.0, 10.0, 10.0])
        state.info.cutoff = 4.5
        
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 创建force
        force = pygcmc.DrudeForce()
        
        charge = -1.71636
        k_spring = 418400.0
        polarizability = 1.71636**2 * 138.935456 / k_spring
        
        for i in range(2):
            force.addParticle(
                drudeIndex=5*i+1, parentIndex=5*i,
                aniso1Index=-1, aniso2Index=-1,
                aniso3Index=-1, aniso4Index=-1,
                charge=charge, polarizability=polarizability,
                aniso12=1.0, aniso34=1.0
            )
        
        force.addScreenedPair(0, 1, 1.3)
        
        # 测试SCF
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.01
        params.maxIterations = 500
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 重置Drude
        for i in range(2):
            state.atoms[5*i+1].x = state.atoms[5*i].x
            state.atoms[5*i+1].y = state.atoms[5*i].y
            state.atoms[5*i+1].z = state.atoms[5*i].z
        
        energy_scf = force.calculateEnergySCF(state)
        
        # 保存SCF位置
        scf_positions = []
        for i in range(2):
            scf_positions.append([
                state.atoms[5*i+1].x,
                state.atoms[5*i+1].y,
                state.atoms[5*i+1].z
            ])
        
        # 测试FBP
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置Drude
        for i in range(2):
            state.atoms[5*i+1].x = state.atoms[5*i].x
            state.atoms[5*i+1].y = state.atoms[5*i].y
            state.atoms[5*i+1].z = state.atoms[5*i].z
        
        energy_fbp = force.calculateEnergySCF(state)
        
        diff = abs(energy_scf - energy_fbp)
        
        print(f"{dist:<10.2f} {energy_scf:<15.6f} {energy_fbp:<15.6f} {diff:<10.6f}")
    
    # 分析能量组成
    print("\n\n能量组成分析（间距0.5nm）：")
    print("-"*60)
    
    # 使用0.5nm的系统
    dist = 0.5
    atoms = []
    residues = []
    
    for i in range(2):
        x_base = i * dist
        positions = [
            (x_base, 0.0, 0.0, 1.71636, 0),   # O
            (x_base, 0.0, 0.0, -1.71636, 1),  # D
            (x_base + 0.09572, 0.0, 0.0, 0.55733, 2),  # H1
            (x_base - 0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
            (x_base, -0.024034, 0.0, -1.11466, 3)  # M
        ]
        
        for x, y, z, charge, typ in positions:
            atom = pygcmc.MCAtom()
            atom.x = x
            atom.y = y
            atom.z = z
            atom.charge = charge
            atom.type = typ
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
    state.activeAtomCount = 10
    state.activeResidueCount = 2
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print("\n计算能量时包含：")
    print("1. Drude谐振子能量")
    print("2. Thole屏蔽的Drude-Drude相互作用")
    print("3. 所有分子间的库仑相互作用")
    print("4. 注意：分子内相互作用应被排除")
    
    # 手动估算能量
    print("\n手动估算（粗略）：")
    # 两个O原子间的相互作用
    r_OO = dist
    E_OO = 138.935 * 1.71636 * 1.71636 / r_OO
    print(f"O-O相互作用: ~{E_OO:.1f} kJ/mol")
    
    # 还有很多其他相互作用...
    print("还有O-H, O-M, H-H, H-M, M-M等相互作用...")
    print("以及Drude的贡献...")

if __name__ == "__main__":
    test_simple_energy()