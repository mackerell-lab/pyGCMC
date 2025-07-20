#!/usr/bin/env python3
"""
诊断FBP算法：为什么力平衡了但能量不对？
"""

import pygcmc
import numpy as np

def create_single_water():
    """创建单个水分子用于测试"""
    atoms = []
    residues = []
    
    # 单个水分子
    positions = [
        (0.0, 0.0, 0.0, 1.71636, 0),   # O
        (0.0, 0.0, 0.0, -1.71636, 1),  # D
        (0.09572, 0.0, 0.0, 0.55733, 2),  # H1
        (-0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
        (0.0, -0.024034, 0.0, -1.11466, 3)  # M
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
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    residues.append(res)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = 1
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def create_two_water():
    """创建两个水分子"""
    atoms = []
    residues = []
    
    # 两个水分子，相距1nm
    for i in range(2):
        base_x = i * 1.0
        positions = [
            (base_x, 0.0, 0.0, 1.71636, 0),   # O
            (base_x, 0.0, 0.0, -1.71636, 1),  # D
            (base_x + 0.09572, 0.0, 0.0, 0.55733, 2),  # H1
            (base_x - 0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
            (base_x, -0.024034, 0.0, -1.11466, 3)  # M
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
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = 2
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def analyze_energy_components(state, force):
    """分析能量的各个组成部分"""
    # 计算总能量
    total_energy = force.calculateEnergySCF(state)
    
    # 获取Drude位置
    drude_positions = []
    parent_positions = []
    for i in range(force.getNumParticles()):
        drude_idx = 5*i + 1
        parent_idx = 5*i
        drude_positions.append([
            state.atoms[drude_idx].x,
            state.atoms[drude_idx].y,
            state.atoms[drude_idx].z
        ])
        parent_positions.append([
            state.atoms[parent_idx].x,
            state.atoms[parent_idx].y,
            state.atoms[parent_idx].z
        ])
    
    # 计算谐振子能量
    harmonic_energy = 0.0
    k = 418400.0  # kJ/mol/nm^2
    for i, (d_pos, p_pos) in enumerate(zip(drude_positions, parent_positions)):
        dx = d_pos[0] - p_pos[0]
        dy = d_pos[1] - p_pos[1]
        dz = d_pos[2] - p_pos[2]
        r2 = dx*dx + dy*dy + dz*dz
        harmonic_energy += 0.5 * k * r2
    
    return total_energy, harmonic_energy, drude_positions

def test_single_water():
    """测试单个水分子"""
    print("测试1：单个水分子（无分子间相互作用）")
    print("="*60)
    
    state = create_single_water()
    
    # Drude参数
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    # 测试不同算法
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        print(f"\n{algo}算法:")
        
        force = pygcmc.DrudeForce()
        force.addParticle(
            drudeIndex=1, parentIndex=0,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001 if algo == "SCF" else 0.1
        params.maxIterations = 500 if algo == "SCF" else 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置Drude位置
        state.atoms[1].x = state.atoms[0].x
        state.atoms[1].y = state.atoms[0].y
        state.atoms[1].z = state.atoms[0].z
        
        # 计算能量
        total_energy, harmonic_energy, drude_pos = analyze_energy_components(state, force)
        
        # Drude位移
        dx = drude_pos[0][0] - state.atoms[0].x
        dy = drude_pos[0][1] - state.atoms[0].y
        dz = drude_pos[0][2] - state.atoms[0].z
        displacement = np.sqrt(dx*dx + dy*dy + dz*dz)
        
        print(f"  总能量: {total_energy:.6f} kJ/mol")
        print(f"  谐振子能量: {harmonic_energy:.6f} kJ/mol")
        print(f"  其他能量: {total_energy - harmonic_energy:.6f} kJ/mol")
        print(f"  Drude位移: {displacement*1000:.3f} pm")

def test_force_vs_energy():
    """测试力和能量的关系"""
    print("\n\n测试2：力-能量关系")
    print("="*60)
    
    state = create_two_water()
    
    # 创建force对象
    force = pygcmc.DrudeForce()
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    for i in range(2):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    force.addScreenedPair(0, 1, 1.3)
    
    # 测试不同Drude位移下的能量
    print("\n手动设置Drude位移，观察能量变化:")
    print(f"{'Displacement (pm)':<20} {'Energy (kJ/mol)':<20} {'dE/dr estimate':<20}")
    print("-"*60)
    
    displacements = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]  # pm
    energies = []
    
    for disp_pm in displacements:
        disp = disp_pm / 1000.0  # pm to nm
        
        # 设置Drude位置（沿x方向位移）
        for i in range(2):
            state.atoms[5*i+1].x = state.atoms[5*i].x + disp
            state.atoms[5*i+1].y = state.atoms[5*i].y
            state.atoms[5*i+1].z = state.atoms[5*i].z
        
        # 计算能量（不优化位置）
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1000.0  # 很大的容差，避免优化
        params.maxIterations = 1
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        energy = force.calculateEnergySCF(state)
        energies.append(energy)
        
        # 估算力（能量梯度）
        if len(energies) > 1:
            de_dr = (energies[-1] - energies[-2]) / (displacements[-1] - displacements[-2])
        else:
            de_dr = 0.0
        
        print(f"{disp_pm:<20.1f} {energy:<20.6f} {de_dr:<20.3f}")

def test_convergence_path():
    """测试收敛路径"""
    print("\n\n测试3：收敛路径分析")
    print("="*60)
    
    state = create_two_water()
    
    # 创建两个force对象用于比较
    force_scf = pygcmc.DrudeForce()
    force_fbp = pygcmc.DrudeForce()
    
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    for force in [force_scf, force_fbp]:
        for i in range(2):
            force.addParticle(
                drudeIndex=5*i+1, parentIndex=5*i,
                aniso1Index=-1, aniso2Index=-1,
                aniso3Index=-1, aniso4Index=-1,
                charge=charge, polarizability=polarizability,
                aniso12=1.0, aniso34=1.0
            )
        force.addScreenedPair(0, 1, 1.3)
    
    # SCF设置
    params_scf = pygcmc.DrudeSCFParams()
    params_scf.tolerance = 0.001
    params_scf.maxIterations = 500
    force_scf.setSCFParameters(params_scf)
    force_scf.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # FBP设置
    params_fbp = pygcmc.DrudeSCFParams()
    params_fbp.tolerance = 0.1
    params_fbp.maxIterations = 50
    force_fbp.setSCFParameters(params_fbp)
    force_fbp.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 重置并运行SCF
    for i in range(2):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    energy_scf = force_scf.calculateEnergySCF(state)
    
    # 保存SCF位置
    scf_positions = []
    for i in range(2):
        scf_positions.append([
            state.atoms[5*i+1].x - state.atoms[5*i].x,
            state.atoms[5*i+1].y - state.atoms[5*i].y,
            state.atoms[5*i+1].z - state.atoms[5*i].z
        ])
    
    # 重置并运行FBP
    for i in range(2):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    energy_fbp = force_fbp.calculateEnergySCF(state)
    
    # 保存FBP位置
    fbp_positions = []
    for i in range(2):
        fbp_positions.append([
            state.atoms[5*i+1].x - state.atoms[5*i].x,
            state.atoms[5*i+1].y - state.atoms[5*i].y,
            state.atoms[5*i+1].z - state.atoms[5*i].z
        ])
    
    print(f"SCF能量: {energy_scf:.6f} kJ/mol")
    print(f"FBP能量: {energy_fbp:.6f} kJ/mol")
    print(f"能量差: {abs(energy_fbp - energy_scf):.6f} kJ/mol")
    
    print("\nDrude位置比较:")
    for i in range(2):
        scf_disp = np.sqrt(sum(x**2 for x in scf_positions[i]))
        fbp_disp = np.sqrt(sum(x**2 for x in fbp_positions[i]))
        print(f"  水分子{i}: SCF位移={scf_disp*1000:.3f}pm, FBP位移={fbp_disp*1000:.3f}pm")

def main():
    """主函数"""
    print("FBP算法诊断：为什么力平衡但能量不对？")
    print("="*80)
    
    test_single_water()
    test_force_vs_energy()
    test_convergence_path()
    
    print("\n\n分析总结:")
    print("="*60)
    print("1. 检查单个水分子的能量分解")
    print("2. 检查力-能量关系是否一致")
    print("3. 比较SCF和FBP的收敛路径")

if __name__ == "__main__":
    main()