#!/usr/bin/env python3
"""
对比OpenMM和PyGCMC的Drude SCF结果
使用相同的水系统，比较：
1. Drude位置
2. 系统能量
3. 诱导偶极矩
"""

import numpy as np
import pickle
import os

try:
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    HAS_OPENMM = True
    print("OpenMM导入成功")
except ImportError:
    HAS_OPENMM = False
    print("警告：OpenMM未安装")

import pygcmc

def create_openmm_system_from_data(data, n_test_waters=10):
    """
    从PyGCMC数据创建OpenMM系统
    """
    n_waters = min(n_test_waters, data['n_waters'])
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"\n创建OpenMM系统 ({n_waters}个水分子)...")
    
    # 创建系统
    system = mm.System()
    
    # 添加粒子
    masses = [15.99943, 0.4, 1.007947, 1.007947, 0.0]  # O, D, H1, H2, M
    
    for i in range(n_waters * 5):
        system.addParticle(masses[i % 5] * unit.amu)
    
    # 设置盒子
    system.setDefaultPeriodicBoxVectors(
        [box_length, 0, 0] * unit.nanometer,
        [0, box_length, 0] * unit.nanometer,
        [0, 0, box_length] * unit.nanometer
    )
    
    # 创建拓扑（简化版）
    topology = app.Topology()
    chain = topology.addChain()
    
    for i in range(n_waters):
        residue = topology.addResidue('HOH', chain)
        o_atom = topology.addAtom('O', app.element.oxygen, residue)
        h1_atom = topology.addAtom('H1', app.element.hydrogen, residue)
        h2_atom = topology.addAtom('H2', app.element.hydrogen, residue)
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)
    
    topology.setPeriodicBoxVectors([
        [box_length, 0, 0],
        [0, box_length, 0],
        [0, 0, box_length]
    ] * unit.nanometer)
    
    # 非键相互作用
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(min(0.9, box_length/2 - 0.01) * unit.nanometer)
    
    # SWM4-NDP参数
    charges = data['charges']
    
    for i in range(n_waters):
        for j in range(5):
            idx = i * 5 + j
            charge = charges[j] * unit.elementary_charge
            sigma = 0.318395 * unit.nanometer if j == 0 else 1.0 * unit.nanometer
            epsilon = 0.88257 * unit.kilojoule_per_mole if j == 0 else 0.0 * unit.kilojoule_per_mole
            
            nonbonded.addParticle(charge, sigma, epsilon)
        
        # 分子内排除
        base = i * 5
        for j in range(5):
            for k in range(j+1, 5):
                nonbonded.addException(base+j, base+k, 0, 1, 0)
    
    system.addForce(nonbonded)
    
    # DrudeForce
    drudeForce = mm.DrudeForce()
    
    for i in range(n_waters):
        parent_idx = i * 5
        drude_idx = i * 5 + 1
        
        drudeForce.addParticle(
            drude_idx,    # Drude
            parent_idx,   # parent
            -1, -1, -1,   # aniso
            -1.71636 * unit.elementary_charge,
            0.0009782237 * unit.nanometer**3,
            1.0, 1.0
        )
    
    # 添加Thole屏蔽
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drudeForce.addScreenedPair(i, j, 1.3)
    
    system.addForce(drudeForce)
    
    # 准备位置（只取前n_waters个水）
    positions_openmm = []
    for i in range(n_waters * 5):
        positions_openmm.append(positions[i])
    positions_openmm = positions_openmm * unit.nanometer
    
    return system, topology, positions_openmm, drudeForce

def run_openmm_scf(system, topology, positions):
    """
    使用OpenMM的SCF积分器优化Drude位置
    """
    print("\n运行OpenMM SCF优化...")
    
    # 创建SCF积分器
    scf_integrator = mm.DrudeSCFIntegrator(0.001 * unit.picosecond)
    # 注意：DrudeSCFIntegrator不支持setMaxDrudeDistance
    
    # 创建模拟
    simulation = app.Simulation(topology, system, scf_integrator)
    
    # 确保所有Drude粒子初始在parent位置
    positions_copy = []
    for i in range(len(positions)):
        if i % 5 == 1:  # Drude粒子
            parent_idx = i - 1
            positions_copy.append(positions[parent_idx])
        else:
            positions_copy.append(positions[i])
    
    simulation.context.setPositions(positions_copy)
    
    # 获取初始能量
    state_initial = simulation.context.getState(getEnergy=True)
    energy_initial = state_initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"  初始能量（Drude在parent位置）: {energy_initial:.2f} kJ/mol")
    
    # 运行SCF（只需要一步）
    simulation.step(1)
    
    # 获取优化后的状态
    state_final = simulation.context.getState(getPositions=True, getEnergy=True)
    positions_final = state_final.getPositions(asNumpy=True)
    energy_final = state_final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    print(f"  优化后能量: {energy_final:.2f} kJ/mol")
    print(f"  能量变化: {energy_final - energy_initial:.2f} kJ/mol")
    
    return positions_final, energy_final

def run_pygcmc_scf(data, n_test_waters=10):
    """
    使用PyGCMC运行SCF优化
    """
    n_waters = min(n_test_waters, data['n_waters'])
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"\n运行PyGCMC SCF优化 ({n_waters}个水分子)...")
    
    # 创建PyGCMC状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = data['charges']
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters * 5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = atom_types[i % 5]
        atoms.append(atom)
    
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(0.9, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    for i in range(n_waters):
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
    
    # 添加Thole对（与OpenMM相同）
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0  # 严格容差
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 确保Drude在parent位置
    state_copy = state.copy()
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        state_copy.atoms[d_idx].x = state_copy.atoms[o_idx].x
        state_copy.atoms[d_idx].y = state_copy.atoms[o_idx].y
        state_copy.atoms[d_idx].z = state_copy.atoms[o_idx].z
    
    # 计算初始能量（粗略估计）
    state_initial = state_copy.copy()
    params_no_opt = pygcmc.DrudeSCFParams()
    params_no_opt.tolerance = 1e10
    params_no_opt.maxIterations = 1
    params_no_opt.dampingFactor = 0.0
    params_no_opt.maxDrudeDistance = 0.00001
    force.setSCFParameters(params_no_opt)
    
    try:
        energy_initial = force.calculateEnergySCF(state_initial)
        print(f"  初始能量（Drude在parent位置）: {energy_initial:.2f} kJ/mol")
    except:
        energy_initial = None
    
    # 运行完整SCF
    force.setSCFParameters(params)
    energy_final = force.calculateEnergySCF(state_copy)
    
    print(f"  优化后能量: {energy_final:.2f} kJ/mol")
    if energy_initial is not None:
        print(f"  能量变化: {energy_final - energy_initial:.2f} kJ/mol")
    
    return state_copy, energy_final

def compare_results(positions_openmm, energy_openmm, state_pygcmc, energy_pygcmc, n_waters):
    """
    对比OpenMM和PyGCMC的结果
    """
    print("\n" + "="*70)
    print("对比OpenMM和PyGCMC的SCF结果")
    print("="*70)
    
    # 能量对比
    print(f"\n能量对比:")
    print(f"  OpenMM: {energy_openmm:.2f} kJ/mol ({energy_openmm/n_waters:.2f} kJ/mol/水)")
    print(f"  PyGCMC: {energy_pygcmc:.2f} kJ/mol ({energy_pygcmc/n_waters:.2f} kJ/mol/水)")
    print(f"  差异: {abs(energy_openmm - energy_pygcmc):.2f} kJ/mol")
    
    # Drude位移对比
    print(f"\nDrude位移对比 (前5个水分子):")
    print(f"{'水':>4} {'OpenMM(pm)':>12} {'PyGCMC(pm)':>12} {'差异(pm)':>10}")
    print("-"*40)
    
    displacements_openmm = []
    displacements_pygcmc = []
    
    for i in range(min(5, n_waters)):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        # OpenMM位移
        o_pos_omm = positions_openmm[o_idx].value_in_unit(unit.nanometer)
        d_pos_omm = positions_openmm[d_idx].value_in_unit(unit.nanometer)
        disp_omm = np.linalg.norm(d_pos_omm - o_pos_omm) * 1000
        displacements_openmm.append(disp_omm)
        
        # PyGCMC位移
        dx = state_pygcmc.atoms[d_idx].x - state_pygcmc.atoms[o_idx].x
        dy = state_pygcmc.atoms[d_idx].y - state_pygcmc.atoms[o_idx].y
        dz = state_pygcmc.atoms[d_idx].z - state_pygcmc.atoms[o_idx].z
        disp_pgc = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        displacements_pygcmc.append(disp_pgc)
        
        diff = abs(disp_omm - disp_pgc)
        
        print(f"{i+1:4d} {disp_omm:12.2f} {disp_pgc:12.2f} {diff:10.2f}")
    
    # 统计分析
    avg_omm = np.mean(displacements_openmm)
    avg_pgc = np.mean(displacements_pygcmc)
    
    print(f"\n平均位移:")
    print(f"  OpenMM: {avg_omm:.2f} pm")
    print(f"  PyGCMC: {avg_pgc:.2f} pm")
    print(f"  差异: {abs(avg_omm - avg_pgc):.2f} pm")
    
    # 诱导偶极矩对比
    print(f"\n诱导偶极矩对比 (前3个水分子):")
    print(f"{'水':>4} {'OpenMM(D)':>10} {'PyGCMC(D)':>10} {'差异(D)':>10}")
    print("-"*35)
    
    q_drude = -1.71636
    
    for i in range(min(3, n_waters)):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        # OpenMM偶极矩
        o_pos_omm = positions_openmm[o_idx].value_in_unit(unit.nanometer)
        d_pos_omm = positions_openmm[d_idx].value_in_unit(unit.nanometer)
        disp_vec_omm = d_pos_omm - o_pos_omm
        dipole_omm = np.linalg.norm(disp_vec_omm * q_drude * 4.80321)
        
        # PyGCMC偶极矩
        dx = state_pygcmc.atoms[d_idx].x - state_pygcmc.atoms[o_idx].x
        dy = state_pygcmc.atoms[d_idx].y - state_pygcmc.atoms[o_idx].y
        dz = state_pygcmc.atoms[d_idx].z - state_pygcmc.atoms[o_idx].z
        disp_vec_pgc = np.array([dx, dy, dz])
        dipole_pgc = np.linalg.norm(disp_vec_pgc * q_drude * 4.80321)
        
        diff = abs(dipole_omm - dipole_pgc)
        
        print(f"{i+1:4d} {dipole_omm:10.3f} {dipole_pgc:10.3f} {diff:10.3f}")
    
    # 结论
    print(f"\n结论:")
    
    energy_diff_percent = abs(energy_openmm - energy_pygcmc) / abs(energy_openmm) * 100
    disp_diff_percent = abs(avg_omm - avg_pgc) / max(avg_omm, avg_pgc) * 100
    
    if energy_diff_percent < 5 and disp_diff_percent < 10:
        print("✓ OpenMM和PyGCMC的SCF结果高度一致")
    elif energy_diff_percent < 10 and disp_diff_percent < 20:
        print("✓ OpenMM和PyGCMC的SCF结果基本一致")
    else:
        print("⚠ OpenMM和PyGCMC的SCF结果存在差异")
        print("  可能的原因：")
        print("  1. 实现细节不同（如收敛判据）")
        print("  2. 数值精度差异")
        print("  3. 算法参数不同")

def main():
    """
    主函数
    """
    print("OpenMM vs PyGCMC Drude SCF对比测试")
    print("="*70)
    
    if not HAS_OPENMM:
        print("错误：需要安装OpenMM")
        return
    
    # 加载水系统
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    # 使用前10个水分子进行测试
    n_test_waters = 10
    
    # 运行OpenMM SCF
    system, topology, positions, drudeForce = create_openmm_system_from_data(data, n_test_waters)
    positions_openmm, energy_openmm = run_openmm_scf(system, topology, positions)
    
    # 运行PyGCMC SCF
    state_pygcmc, energy_pygcmc = run_pygcmc_scf(data, n_test_waters)
    
    # 对比结果
    compare_results(positions_openmm, energy_openmm, state_pygcmc, energy_pygcmc, n_test_waters)

if __name__ == "__main__":
    main()