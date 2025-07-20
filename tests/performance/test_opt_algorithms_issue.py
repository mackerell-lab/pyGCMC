#!/usr/bin/env python3
"""
诊断OPT系列算法的问题
"""

import pygcmc
import numpy as np

def diagnose_opt_algorithms():
    """诊断OPT算法为什么表现这么差"""
    print("OPT算法问题诊断")
    print("="*60)
    
    # 创建简单的2水分子系统
    atoms = []
    residues = []
    
    # 两个水分子
    for i in range(2):
        x_base = i * 0.5
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
    
    state.info.box = np.array([5.0, 5.0, 5.0])
    state.info.cutoff = 2.5
    
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
    
    # 测试各算法
    algorithms = ["SCF", "OPT3", "OPT4", "HybridOPT", "FBP"]
    
    print("\n算法对比：")
    print(f"{'算法':<12} {'能量(kJ/mol)':<15} {'RMSD(pm)':<12} {'诊断':<30}")
    print("-"*70)
    
    for algo in algorithms:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1.0
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        elif algo == "OPT3":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        elif algo == "OPT4":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT4)
        elif algo == "HybridOPT":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.HybridOPT)
        elif algo == "FBP":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置Drude位置
        for i in range(2):
            state.atoms[5*i+1].x = state.atoms[5*i].x
            state.atoms[5*i+1].y = state.atoms[5*i].y
            state.atoms[5*i+1].z = state.atoms[5*i].z
        
        try:
            energy = force.calculateEnergySCF(state)
            
            # 计算RMSD
            displacements = []
            for i in range(2):
                dx = state.atoms[5*i+1].x - state.atoms[5*i].x
                dy = state.atoms[5*i+1].y - state.atoms[5*i].y
                dz = state.atoms[5*i+1].z - state.atoms[5*i].z
                disp = np.sqrt(dx*dx + dy*dy + dz*dz)
                displacements.append(disp)
            
            rmsd = np.sqrt(np.mean(np.array(displacements)**2)) * 1000  # pm
            
            # 诊断
            if energy > 0:
                diagnosis = "能量为正（物理错误）"
            elif rmsd > 5:
                diagnosis = "位移过大"
            elif abs(energy) < 1.0:
                diagnosis = "能量过小（可能未收敛）"
            else:
                diagnosis = "正常"
            
            print(f"{algo:<12} {energy:<15.6f} {rmsd:<12.3f} {diagnosis:<30}")
            
        except Exception as e:
            print(f"{algo:<12} {'ERROR':<15} {'N/A':<12} {str(e)[:30]:<30}")
    
    print("\n\n问题分析：")
    print("-"*60)
    print("OPT3/OPT4/HybridOPT的问题可能是：")
    print("1. 历史位置的权重系数不合适")
    print("2. 没有正确处理初始迭代（历史不足时）")
    print("3. 可能在计算能量时使用了错误的Drude位置")
    print("4. HybridOPT可能没有正确执行最后的SCF步骤")
    
    # 测试OPT3的收敛过程
    print("\n\n测试OPT3的迭代过程：")
    print("-"*60)
    
    # 尝试手动设置一些中间位置
    force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    
    # 给Drude一些初始位移
    for i in range(2):
        state.atoms[5*i+1].x = state.atoms[5*i].x + 0.001  # 1 pm
    
    print("初始位移后的能量：")
    energy = force.calculateEnergySCF(state)
    print(f"能量: {energy:.6f} kJ/mol")
    
    # 检查最终位移
    final_disp = []
    for i in range(2):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        final_disp.append(dx*1000)  # pm
    print(f"最终位移: {final_disp}")

if __name__ == "__main__":
    diagnose_opt_algorithms()