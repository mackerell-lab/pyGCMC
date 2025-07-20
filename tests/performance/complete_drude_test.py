#!/usr/bin/env python3
"""
完整的Drude测试，包含所有必要的能量项
"""

import numpy as np
import pygcmc

try:
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

def test_pygcmc_with_full_energy():
    """
    测试PyGCMC的Drude+完整能量
    """
    print("\nPyGCMC完整能量测试")
    print("="*70)
    
    # 创建简单的两水系统
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
        [0.3, 0.3, 0.3],      # D (初始在O位置)
        [0.396, 0.3, 0.3],    # H1
        [0.252, 0.377, 0.3],  # H2
        [0.3, 0.3, 0.3],      # M
        # 水2
        [0.7, 0.7, 0.7],      # O
        [0.7, 0.7, 0.7],      # D (初始在O位置)
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
    
    # 创建4x4的参数矩阵（只有O-O相互作用有非零值）
    ljSigma = []
    ljEps = []
    for i in range(4):
        for j in range(4):
            if i == 0 and j == 0:  # O-O相互作用
                ljSigma.append(0.318395)
                ljEps.append(0.88257)
            else:
                ljSigma.append(0.0)
                ljEps.append(0.0)
    
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    # 1. 先计算非Drude能量（库仑+LJ）
    print("\n1. 非Drude能量计算")
    
    # 设置PME参数
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]  # 适合小系统的网格
    spline_order = 4
    
    # 初始化PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(
        state.info.cutoff, 
        [state.info.box[0], state.info.box[1], state.info.box[2]], 
        alpha, 
        mesh_size, 
        spline_order
    )
    
    # 计算能量
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    energy_no_drude = pme_elec + pme_vdw
    print(f"  库仑+LJ能量: {energy_no_drude:.4f} kJ/mol")
    print(f"    库仑: {pme_elec:.4f} kJ/mol")
    print(f"    LJ: {pme_vdw:.4f} kJ/mol")
    
    # 2. Drude SCF优化
    print("\n2. Drude SCF优化")
    
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
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 运行SCF
    state_opt = state.copy()
    
    try:
        # 注意：calculateEnergySCF只返回Drude相关的能量
        energy_drude = force.calculateEnergySCF(state_opt)
        print(f"  Drude能量: {energy_drude:.4f} kJ/mol")
        
        # 分析Drude位移
        print("\n  Drude位移分析:")
        for i in range(2):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state_opt.atoms[d_idx].x - state_opt.atoms[o_idx].x
            dy = state_opt.atoms[d_idx].y - state_opt.atoms[o_idx].y
            dz = state_opt.atoms[d_idx].z - state_opt.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            
            # 计算诱导偶极矩
            q_drude = -1.71636
            dipole_vec = np.array([dx, dy, dz]) * q_drude * 4.80321  # Debye
            dipole_mag = np.linalg.norm(dipole_vec)
            
            print(f"    水{i+1}: 位移 = {disp:.2f} pm, 偶极矩 = {dipole_mag:.3f} D")
        
        # 3. 计算优化后的总能量
        print("\n3. 优化后的总能量")
        pme_elec_opt, pme_vdw_opt, pme_dict_opt = pygcmc.computeSystemEnergyPME(state_opt)
        energy_total_opt = pme_elec_opt + pme_vdw_opt
        print(f"  优化后库仑+LJ: {energy_total_opt:.4f} kJ/mol")
        print(f"    库仑: {pme_elec_opt:.4f} kJ/mol")
        print(f"    LJ: {pme_vdw_opt:.4f} kJ/mol")
        print(f"  总能量变化: {energy_total_opt - energy_no_drude:.4f} kJ/mol")
        
    except Exception as e:
        print(f"  SCF失败: {e}")

def test_openmm_with_full_energy():
    """
    测试OpenMM的完整Drude系统
    """
    if not HAS_OPENMM:
        return
        
    print("\n\nOpenMM完整能量测试")
    print("="*70)
    
    # 创建系统
    system = mm.System()
    
    # 参数
    box_size = 1.0  # nm
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    masses = [15.99943, 0.4, 1.007947, 1.007947, 0.0]
    
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
    
    # 添加粒子
    for i in range(10):
        system.addParticle(masses[i % 5] * unit.amu)
    
    # 设置盒子
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0] * unit.nanometer,
        [0, box_size, 0] * unit.nanometer,
        [0, 0, box_size] * unit.nanometer
    )
    
    # 1. NonbondedForce
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(0.45 * unit.nanometer)
    
    for i in range(10):
        charge = charges[i % 5] * unit.elementary_charge
        sigma = 0.318395 * unit.nanometer if i % 5 == 0 else 1.0 * unit.nanometer
        epsilon = 0.88257 * unit.kilojoule_per_mole if i % 5 == 0 else 0.0 * unit.kilojoule_per_mole
        nonbonded.addParticle(charge, sigma, epsilon)
    
    # 分子内排除
    for i in range(2):
        base = i * 5
        for j in range(5):
            for k in range(j+1, 5):
                nonbonded.addException(base+j, base+k, 0, 1, 0)
    
    system.addForce(nonbonded)
    
    # 2. DrudeForce
    drudeForce = mm.DrudeForce()
    
    for i in range(2):
        drudeForce.addParticle(
            i*5+1,    # Drude
            i*5,      # parent
            -1, -1, -1,
            -1.71636 * unit.elementary_charge,
            0.0009782237 * unit.nanometer**3,
            1.0, 1.0
        )
    
    drudeForce.addScreenedPair(0, 1, 1.3)
    system.addForce(drudeForce)
    
    # 初始能量（Drude在parent位置）
    print("\n1. 初始能量（Drude在parent位置）")
    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(system, integrator)
    context.setPositions(positions * unit.nanometer)
    
    state = context.getState(getEnergy=True)
    energy_init = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"  初始能量: {energy_init:.4f} kJ/mol")
    
    # SCF优化
    print("\n2. SCF优化")
    scf_integrator = mm.DrudeSCFIntegrator(0.001 * unit.picosecond)
    scf_context = mm.Context(system, scf_integrator)
    scf_context.setPositions(positions * unit.nanometer)
    
    # 运行SCF
    scf_integrator.step(1)
    
    state_scf = scf_context.getState(getPositions=True, getEnergy=True)
    energy_scf = state_scf.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    positions_scf = state_scf.getPositions(asNumpy=True)
    
    print(f"  SCF后能量: {energy_scf:.4f} kJ/mol")
    print(f"  能量变化: {energy_scf - energy_init:.4f} kJ/mol")
    
    # 分析Drude位移
    print("\n  Drude位移分析:")
    for i in range(2):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        o_pos = positions_scf[o_idx].value_in_unit(unit.nanometer)
        d_pos = positions_scf[d_idx].value_in_unit(unit.nanometer)
        
        disp = np.linalg.norm(d_pos - o_pos) * 1000
        
        # 计算诱导偶极矩
        q_drude = -1.71636
        dipole_vec = (d_pos - o_pos) * q_drude * 4.80321  # Debye
        dipole_mag = np.linalg.norm(dipole_vec)
        
        print(f"    水{i+1}: 位移 = {disp:.2f} pm, 偶极矩 = {dipole_mag:.3f} D")

def main():
    """
    主函数
    """
    print("完整Drude系统测试")
    print("="*70)
    
    # 测试PyGCMC
    test_pygcmc_with_full_energy()
    
    # 测试OpenMM
    test_openmm_with_full_energy()
    
    print("\n\n结论:")
    print("="*70)
    print("1. Drude粒子需要外电场才能产生位移")
    print("2. 电场来自其他带电粒子（通过NonbondedForce）")
    print("3. PyGCMC和OpenMM的能量计算方式不同")
    print("4. 需要包含所有能量项才能正确比较")

if __name__ == "__main__":
    main()