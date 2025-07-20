#!/usr/bin/env python3
"""
测试PyGCMC的完整Drude能量计算
使用computeSystemEnergyDrude
"""

import numpy as np
import pygcmc

def test_drude_complete_energy():
    """
    测试完整的Drude能量计算
    """
    print("PyGCMC完整Drude能量测试")
    print("="*70)
    
    # 创建两个水分子系统
    state = pygcmc.MCState()
    
    box_size = 1.0  # nm
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    # 位置
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
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    
    ljSigma = []
    ljEps = []
    for i in range(4):
        for j in range(4):
            if i == 0 and j == 0:
                ljSigma.append(0.318395)
                ljEps.append(0.88257)
            else:
                ljSigma.append(0.0)
                ljEps.append(0.0)
    
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    # 初始化Drude系统
    print("\n1. 初始化Drude系统")
    
    # 清除之前的设置
    pygcmc.clearDrudeForce()
    
    # 先初始化
    pygcmc.initializeDrudeForce()
    
    # 添加Drude粒子（使用全局函数）
    for i in range(2):
        pygcmc.addDrudeParticle(
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
    
    # 添加Thole屏蔽
    pygcmc.addDrudeScreenedPair(0, 1, 1.3)
    
    print(f"  添加了 {pygcmc.getNumDrudeParticles()} 个Drude粒子")
    print(f"  添加了 {pygcmc.getNumDrudeScreenedPairs()} 个Thole对")
    
    # 设置SCF参数
    pygcmc.setDrudeSCFTolerance(1.0)  # kJ/mol
    
    # 2. 计算初始能量（Drude在parent位置）
    print("\n2. 初始能量（Drude在parent位置）")
    
    try:
        energy_result = pygcmc.computeSystemEnergyDrude(state)
        print(f"  结果类型: {type(energy_result)}")
        
        if isinstance(energy_result, tuple):
            if len(energy_result) >= 3:
                energy_elec, energy_vdw, energy_dict = energy_result
                print(f"  总库仑能量: {energy_elec:.4f} kJ/mol")
                print(f"  总LJ能量: {energy_vdw:.4f} kJ/mol")
                print(f"  总能量: {energy_elec + energy_vdw:.4f} kJ/mol")
                
                if isinstance(energy_dict, dict):
                    print("\n  能量组成:")
                    for key, value in energy_dict.items():
                        print(f"    {key}: {value:.4f} kJ/mol")
            else:
                print(f"  返回值: {energy_result}")
        else:
            print(f"  能量: {energy_result:.4f} kJ/mol")
            
    except Exception as e:
        print(f"  计算失败: {e}")
        import traceback
        traceback.print_exc()
    
    # 3. 检查Drude位置
    print("\n3. 检查优化后的Drude位置")
    
    # 看看state是否被修改了
    for i in range(2):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = state.atoms[d_idx].x - state.atoms[o_idx].x
        dy = state.atoms[d_idx].y - state.atoms[o_idx].y
        dz = state.atoms[d_idx].z - state.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        
        print(f"  水{i+1} Drude位移: {disp:.4f} pm")
    
    # 4. 手动测试DrudeForce
    print("\n4. 使用DrudeForce类测试")
    
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
    
    # 设置使用computeSystemEnergyDrude
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 测试能量计算
    state_copy = state.copy()
    energy_drude_only = force.calculateEnergySCF(state_copy)
    print(f"  DrudeForce能量: {energy_drude_only:.4f} kJ/mol")
    
    # 检查位移
    for i in range(2):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = state_copy.atoms[d_idx].x - state_copy.atoms[o_idx].x
        dy = state_copy.atoms[d_idx].y - state_copy.atoms[o_idx].y
        dz = state_copy.atoms[d_idx].z - state_copy.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        
        print(f"  水{i+1} Drude位移: {disp:.4f} pm")

def main():
    """
    主函数
    """
    test_drude_complete_energy()
    
    print("\n\n总结:")
    print("="*70)
    print("computeSystemEnergyDrude应该计算完整的带Drude的系统能量")
    print("包括：库仑相互作用 + LJ相互作用 + Drude极化")

if __name__ == "__main__":
    main()