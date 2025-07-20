#!/usr/bin/env python3
"""
分析为什么CG给出更低的能量
"""

import pygcmc
import numpy as np

def create_simple_system(n_waters=8):
    """创建简单的水系统用于分析"""
    spacing = 0.5
    atoms = []
    residues = []
    
    for i in range(n_waters):
        x = (i % 2) * spacing + 0.5
        y = (i // 2) * spacing + 0.5
        z = 0.5
        
        # SWM4-NDP水模型
        positions = [
            (x, y, z, 1.71636, 0),   # O
            (x, y, z, -1.71636, 1),  # D
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
    
    box_length = 3.0
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = 1.2
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def calculate_detailed_energy(state, force, n_waters):
    """计算详细的能量分解"""
    # 计算各部分能量
    
    # 1. Drude弹簧能量
    spring_energy = 0.0
    for i in range(n_waters):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        dy = state.atoms[5*i+1].y - state.atoms[5*i].y
        dz = state.atoms[5*i+1].z - state.atoms[5*i].z
        r2 = dx*dx + dy*dy + dz*dz
        
        # E = 0.5 * k * r^2
        k_spring = 418400.0  # kJ/mol/nm^2
        spring_energy += 0.5 * k_spring * r2
    
    # 2. 计算库仑能量（简化）
    coulomb_energy = 0.0
    COULOMB_CONSTANT = 138.935456
    
    for i in range(state.activeAtomCount):
        for j in range(i+1, state.activeAtomCount):
            # 检查是否同一分子
            res_i = -1
            res_j = -1
            for r in range(n_waters):
                if i >= 5*r and i < 5*(r+1):
                    res_i = r
                if j >= 5*r and j < 5*(r+1):
                    res_j = r
            
            if res_i == res_j:
                continue  # 跳过分子内相互作用
            
            dx = state.atoms[i].x - state.atoms[j].x
            dy = state.atoms[i].y - state.atoms[j].y
            dz = state.atoms[i].z - state.atoms[j].z
            
            # PBC
            if state.info.box[0] > 0:
                dx -= state.info.box[0] * round(dx / state.info.box[0])
                dy -= state.info.box[1] * round(dy / state.info.box[1])
                dz -= state.info.box[2] * round(dz / state.info.box[2])
            
            r = np.sqrt(dx*dx + dy*dy + dz*dz)
            if r > 1e-6 and r < state.info.cutoff:
                coulomb_energy += COULOMB_CONSTANT * state.atoms[i].charge * state.atoms[j].charge / r
    
    return {
        'spring': spring_energy,
        'coulomb': coulomb_energy,
        'total': spring_energy + coulomb_energy
    }

def analyze_convergence():
    """分析SCF和CG的收敛行为"""
    print("分析SCF和CG的能量最小化")
    print("="*60)
    
    n_waters = 8
    state = create_simple_system(n_waters)
    
    # 创建力对象
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
    
    # 添加所有Thole对
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 测试不同容差
    tolerances = [0.1, 1.0, 10.0, 100.0]
    
    print(f"\n{'容差':<15} {'算法':<10} {'总能量':<15} {'弹簧能':<15} {'库仑能':<15} {'平均位移(pm)':<15}")
    print("-"*90)
    
    for tol in tolerances:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 500  # 充分迭代
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        # SCF
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        state_scf = state.copy()
        energy_scf = force.calculateEnergySCF(state_scf)
        detailed_scf = calculate_detailed_energy(state_scf, force, n_waters)
        
        # 计算平均位移
        disp_scf = []
        for i in range(n_waters):
            dx = state_scf.atoms[5*i+1].x - state_scf.atoms[5*i].x
            dy = state_scf.atoms[5*i+1].y - state_scf.atoms[5*i].y
            dz = state_scf.atoms[5*i+1].z - state_scf.atoms[5*i].z
            disp_scf.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
        
        print(f"{tol:<15.1f} {'SCF':<10} {energy_scf:<15.2f} {detailed_scf['spring']:<15.2f} "
              f"{detailed_scf['coulomb']:<15.2f} {np.mean(disp_scf):<15.3f}")
        
        # CG
        force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
        state_cg = state.copy()
        energy_cg = force.calculateEnergySCF(state_cg)
        detailed_cg = calculate_detailed_energy(state_cg, force, n_waters)
        
        disp_cg = []
        for i in range(n_waters):
            dx = state_cg.atoms[5*i+1].x - state_cg.atoms[5*i].x
            dy = state_cg.atoms[5*i+1].y - state_cg.atoms[5*i].y
            dz = state_cg.atoms[5*i+1].z - state_cg.atoms[5*i].z
            disp_cg.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
        
        print(f"{'':<15} {'CG':<10} {energy_cg:<15.2f} {detailed_cg['spring']:<15.2f} "
              f"{detailed_cg['coulomb']:<15.2f} {np.mean(disp_cg):<15.3f}")
        
        # 差异
        print(f"{'':<15} {'差异':<10} {energy_cg-energy_scf:<15.2f} "
              f"{detailed_cg['spring']-detailed_scf['spring']:<15.2f} "
              f"{detailed_cg['coulomb']-detailed_scf['coulomb']:<15.2f}")
        print()
    
    # 分析能量最小化原理
    print("\n能量最小化原理分析:")
    print("-"*60)
    print("1. Drude模型的能量函数:")
    print("   E = E_spring + E_coulomb")
    print("   E_spring = Σ 0.5 * k * |r_drude - r_parent|²")
    print("   E_coulomb = Σ q_i * q_j / r_ij")
    print()
    print("2. 最小化目标:")
    print("   找到Drude位置使总能量最小")
    print("   平衡条件: F_spring + F_electric = 0")
    print()
    print("3. 为什么CG能量更低?")
    print("   - CG直接求解线性方程组，保证找到真正的最小值")
    print("   - SCF是迭代方法，可能过早停止或陷入局部最优")
    print("   - 大系统中SCF的收敛性差，导致未充分优化")
    print()
    print("4. 物理含义:")
    print("   - 更低的能量 = 更稳定的构型")
    print("   - CG找到的是更接近真实平衡态的解")
    print("   - 能量差异反映了SCF的收敛不充分")

def test_energy_landscape():
    """测试能量景观"""
    print("\n\n测试单个Drude的能量景观")
    print("="*60)
    
    # 创建只有一个水分子的系统
    state = create_simple_system(1)
    force = pygcmc.DrudeForce()
    
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # 在没有外场的情况下，最优位置应该是r=0
    print("\n1. 无外场情况:")
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    force.setSCFParameters(params)
    
    # 测试不同初始位移
    displacements = [0.0, 0.001, 0.002, 0.005, 0.01]  # nm
    print(f"{'初始位移(pm)':<15} {'弹簧能(kJ/mol)':<20} {'总能量(kJ/mol)':<20}")
    print("-"*55)
    
    for disp in displacements:
        state_test = state.copy()
        state_test.atoms[1].x = state_test.atoms[0].x + disp
        
        # 计算能量
        spring_e = 0.5 * k_spring * disp * disp
        total_e = force.calculateEnergySCF(state_test)
        
        print(f"{disp*1000:<15.1f} {spring_e:<20.6f} {total_e:<20.6f}")
    
    print("\n结论: 能量随位移平方增长，最小值在位移=0处")

if __name__ == "__main__":
    analyze_convergence()
    test_energy_landscape()