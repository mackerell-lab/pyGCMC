#!/usr/bin/env python3
"""
详细对比OpenMM和PyGCMC的Drude实现
分析为什么能量和位移差异巨大
"""

import numpy as np
import pickle
import os

try:
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    print("警告：OpenMM未安装")

import pygcmc

def analyze_drude_parameters():
    """
    分析Drude参数设置
    """
    print("\nDrude参数分析")
    print("="*70)
    
    # SWM4-NDP参数
    print("\nSWM4-NDP水模型参数：")
    print(f"  Drude电荷: -1.71636 e")
    print(f"  极化率: 0.0009782237 nm³")
    
    # 计算力常数
    # k = q²/(2α)
    q_drude = -1.71636  # e
    alpha = 0.0009782237  # nm³
    
    # OpenMM使用的单位转换
    # 1 e² = 138.935456 kJ·mol⁻¹·nm
    ONE_4PI_EPS0 = 138.935456
    
    k_calc = q_drude * q_drude * ONE_4PI_EPS0 / (2 * alpha)
    print(f"\n计算的力常数 k = q²/(2α):")
    print(f"  k = {k_calc:.1f} kJ/(mol·nm²)")
    print(f"  k = {k_calc/100:.1f} kJ/(mol·Å²)")
    
    # OpenMM默认使用的力常数
    k_openmm_default = 418400.0  # kJ/(mol·nm²)
    print(f"\nOpenMM常用力常数:")
    print(f"  k = {k_openmm_default:.1f} kJ/(mol·nm²)")
    
    print(f"\n比值: {k_openmm_default/k_calc:.2f}")

def test_simple_two_water():
    """
    测试简单的两水系统
    """
    print("\n\n简单两水系统测试")
    print("="*70)
    
    # 创建相同的初始位置
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
    
    box_size = 1.0  # nm
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    # 1. PyGCMC测试
    print("\n1. PyGCMC系统")
    print("-"*60)
    
    state = pygcmc.MCState()
    
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
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
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
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
    
    # 计算能量
    state_opt = state.copy()
    try:
        energy_pygcmc = force.calculateEnergySCF(state_opt)
        print(f"  能量: {energy_pygcmc:.4f} kJ/mol")
        
        # 分析Drude位移
        for i in range(2):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state_opt.atoms[d_idx].x - state_opt.atoms[o_idx].x
            dy = state_opt.atoms[d_idx].y - state_opt.atoms[o_idx].y
            dz = state_opt.atoms[d_idx].z - state_opt.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            print(f"  水{i+1} Drude位移: {disp:.2f} pm")
            
    except Exception as e:
        print(f"  计算失败: {e}")
    
    # 2. OpenMM测试
    if HAS_OPENMM:
        print("\n2. OpenMM系统")
        print("-"*60)
        
        system = mm.System()
        
        # 添加粒子
        masses = [15.99943, 0.4, 1.007947, 1.007947, 0.0]
        for i in range(10):
            system.addParticle(masses[i % 5] * unit.amu)
        
        # 设置盒子
        system.setDefaultPeriodicBoxVectors(
            [box_size, 0, 0] * unit.nanometer,
            [0, box_size, 0] * unit.nanometer,
            [0, 0, box_size] * unit.nanometer
        )
        
        # DrudeForce - 使用与PyGCMC相同的参数
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
        
        # 添加NonbondedForce以产生电场
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
        
        # 创建SCF积分器
        scf_integrator = mm.DrudeSCFIntegrator(0.001 * unit.picosecond)
        context = mm.Context(system, scf_integrator)
        context.setPositions(positions * unit.nanometer)
        
        # 运行SCF
        scf_integrator.step(1)
        
        state_omm = context.getState(getPositions=True, getEnergy=True)
        energy_openmm = state_omm.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        positions_omm = state_omm.getPositions(asNumpy=True)
        
        print(f"  能量: {energy_openmm:.4f} kJ/mol")
        
        # 分析Drude位移
        for i in range(2):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            o_pos = positions_omm[o_idx].value_in_unit(unit.nanometer)
            d_pos = positions_omm[d_idx].value_in_unit(unit.nanometer)
            
            disp = np.linalg.norm(d_pos - o_pos) * 1000
            print(f"  水{i+1} Drude位移: {disp:.2f} pm")
    
    # 3. 手动计算预期位移
    print("\n3. 理论分析")
    print("-"*60)
    
    # 计算水1的O原子处的电场
    ONE_4PI_EPS0 = 138.935456
    
    o1_pos = np.array(positions[0])
    field = np.zeros(3)
    
    # 来自水2的电场贡献
    for j in [5, 7, 8, 9]:  # O2, H1, H2, M
        other_pos = np.array(positions[j])
        q_other = charges[j % 5]
        
        delta = o1_pos - other_pos
        r2 = np.dot(delta, delta)
        r = np.sqrt(r2)
        
        E_mag = ONE_4PI_EPS0 * q_other / r2
        field += E_mag * delta / r
    
    field_mag = np.linalg.norm(field)
    print(f"  水1 O原子处的电场: |E| = {field_mag:.1f} kJ/(mol·nm·e)")
    
    # 预期的Drude位移
    alpha = 0.0009782237  # nm³
    q_drude = -1.71636
    expected_disp = alpha * field_mag / abs(q_drude)  # nm
    print(f"  预期Drude位移: {expected_disp*1000:.2f} pm")
    
    # 力常数分析
    k_from_params = q_drude**2 * ONE_4PI_EPS0 / (2 * alpha)
    print(f"\n  从参数计算的力常数: {k_from_params:.1f} kJ/(mol·nm²)")
    
    # 如果使用不同的力常数
    k_openmm = 418400.0
    disp_with_k_openmm = field_mag / k_openmm  # nm
    print(f"  如果k={k_openmm:.0f}: 位移 = {disp_with_k_openmm*1000:.2f} pm")

def analyze_implementation_differences():
    """
    分析实现差异
    """
    print("\n\n实现差异分析")
    print("="*70)
    
    print("\n可能的差异点：")
    print("1. 力常数定义:")
    print("   - PyGCMC可能使用: k = q²/(2α)")
    print("   - OpenMM可能使用固定值或不同公式")
    
    print("\n2. 能量计算:")
    print("   - PyGCMC: 可能包含所有能量项")
    print("   - OpenMM DrudeSCFIntegrator: 可能只计算Drude相关能量")
    
    print("\n3. SCF算法:")
    print("   - 收敛判据不同")
    print("   - 迭代方法不同")
    
    print("\n4. Thole屏蔽:")
    print("   - 实现细节可能不同")
    print("   - 屏蔽函数形式可能不同")

def main():
    """
    主函数
    """
    print("OpenMM和PyGCMC Drude实现详细对比")
    print("="*70)
    
    # 分析参数
    analyze_drude_parameters()
    
    # 测试简单系统
    test_simple_two_water()
    
    # 分析差异
    analyze_implementation_differences()
    
    print("\n\n总结:")
    print("="*70)
    print("1. OpenMM和PyGCMC在Drude实现上存在根本差异")
    print("2. 能量差异可能来自不同的能量项计算")
    print("3. Drude位移差异可能来自不同的力常数定义")
    print("4. 需要查看源代码才能完全理解差异")

if __name__ == "__main__":
    main()