#!/usr/bin/env python3
"""
探索OpenMM的Drude粒子输出功能
"""

import numpy as np

try:
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    print("OpenMM导入成功")
    print(f"OpenMM版本: {mm.__version__}")
except ImportError:
    print("错误：需要安装OpenMM")
    exit(1)

def create_simple_drude_system():
    """
    创建一个简单的Drude系统来探索功能
    """
    print("\n创建简单的2水Drude系统...")
    
    # 创建系统
    system = mm.System()
    
    # SWM4-NDP参数
    mass_O = 15.99943
    mass_H = 1.007947
    mass_D = 0.4  # Drude粒子质量
    mass_M = 0.0  # 虚拟位点
    
    # 添加2个水分子的粒子
    for i in range(2):
        system.addParticle(mass_O * unit.amu)  # O
        system.addParticle(mass_D * unit.amu)  # D
        system.addParticle(mass_H * unit.amu)  # H1
        system.addParticle(mass_H * unit.amu)  # H2
        system.addParticle(mass_M * unit.amu)  # M
    
    # 创建拓扑
    topology = app.Topology()
    chain = topology.addChain()
    
    for i in range(2):
        residue = topology.addResidue('HOH', chain)
        o_atom = topology.addAtom('O', app.element.oxygen, residue)
        # 注意：Drude粒子通常不在拓扑中显示
        h1_atom = topology.addAtom('H1', app.element.hydrogen, residue)
        h2_atom = topology.addAtom('H2', app.element.hydrogen, residue)
        
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)
    
    # 设置盒子
    box_length = 2.0  # nm
    system.setDefaultPeriodicBoxVectors(
        [box_length, 0, 0] * unit.nanometer,
        [0, box_length, 0] * unit.nanometer,
        [0, 0, box_length] * unit.nanometer
    )
    topology.setPeriodicBoxVectors([
        [box_length, 0, 0],
        [0, box_length, 0],
        [0, 0, box_length]
    ] * unit.nanometer)
    
    # 添加非键相互作用
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(0.9 * unit.nanometer)
    
    # SWM4-NDP电荷
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    for i in range(2):  # 2个水分子
        for j in range(5):  # 每个水5个位点
            charge = charges[j] * unit.elementary_charge
            sigma = 0.318395 * unit.nanometer if j == 0 else 1.0 * unit.nanometer
            epsilon = 0.88257 * unit.kilojoule_per_mole if j == 0 else 0.0 * unit.kilojoule_per_mole
            
            nonbonded.addParticle(charge, sigma, epsilon)
        
        # 添加分子内排除
        base = i * 5
        for j in range(5):
            for k in range(j+1, 5):
                nonbonded.addException(base+j, base+k, 0, 1, 0)
    
    system.addForce(nonbonded)
    
    # 添加DrudeForce
    drudeForce = mm.DrudeForce()
    
    k_drude = 418400.0 * unit.kilojoule_per_mole / unit.nanometer**2
    drude_charge = -1.71636 * unit.elementary_charge
    polarizability = 0.00097825258 * unit.nanometer**3
    
    for i in range(2):
        parent_idx = i * 5      # O
        drude_idx = i * 5 + 1   # D
        
        drudeForce.addParticle(
            drude_idx,    # Drude particle
            parent_idx,   # parent particle
            -1, -1, -1,   # aniso indices
            drude_charge,
            polarizability,
            1.0, 1.0      # aniso scales
        )
    
    # 添加Thole屏蔽
    drudeForce.addScreenedPair(0, 1, 1.3)
    
    system.addForce(drudeForce)
    
    return system, topology, drudeForce

def test_drude_positions():
    """
    测试获取Drude粒子位置
    """
    print("\n" + "="*70)
    print("测试Drude粒子位置获取")
    print("="*70)
    
    # 创建系统
    system, topology, drudeForce = create_simple_drude_system()
    
    # 创建初始位置
    positions = []
    # 水1
    positions.extend([
        [0.5, 0.5, 0.5],      # O
        [0.5, 0.5, 0.5],      # D (初始在O位置)
        [0.596, 0.5, 0.5],    # H1
        [0.452, 0.577, 0.5],  # H2
        [0.5, 0.5, 0.5]       # M
    ])
    # 水2
    positions.extend([
        [1.0, 1.0, 1.0],      # O
        [1.0, 1.0, 1.0],      # D (初始在O位置)
        [1.096, 1.0, 1.0],    # H1
        [0.952, 1.077, 1.0],  # H2
        [1.0, 1.0, 1.0]       # M
    ])
    positions = positions * unit.nanometer
    
    # 创建积分器
    print("\n1. 使用DrudeLangevinIntegrator")
    integrator = mm.DrudeLangevinIntegrator(
        300*unit.kelvin,
        1/unit.picosecond,
        1*unit.kelvin,
        10/unit.picosecond,
        0.5*unit.femtosecond
    )
    integrator.setMaxDrudeDistance(0.02 * unit.nanometer)
    
    # 创建模拟
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    
    # 获取初始状态
    print("\n初始Drude位置:")
    state = simulation.context.getState(getPositions=True)
    positions_initial = state.getPositions(asNumpy=True)
    
    for i in range(2):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        o_pos = positions_initial[o_idx].value_in_unit(unit.nanometer)
        d_pos = positions_initial[d_idx].value_in_unit(unit.nanometer)
        
        disp = d_pos - o_pos
        disp_mag = np.linalg.norm(disp) * 1000  # pm
        
        print(f"  水{i+1}: O=({o_pos[0]:.3f}, {o_pos[1]:.3f}, {o_pos[2]:.3f})")
        print(f"       D=({d_pos[0]:.3f}, {d_pos[1]:.3f}, {d_pos[2]:.3f})")
        print(f"       位移={disp_mag:.2f} pm")
    
    # 运行几步模拟
    print("\n运行10步模拟...")
    simulation.step(10)
    
    # 获取新状态
    print("\n模拟后Drude位置:")
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    positions_final = state.getPositions(asNumpy=True)
    energy = state.getPotentialEnergy()
    
    for i in range(2):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        o_pos = positions_final[o_idx].value_in_unit(unit.nanometer)
        d_pos = positions_final[d_idx].value_in_unit(unit.nanometer)
        
        disp = d_pos - o_pos
        disp_mag = np.linalg.norm(disp) * 1000  # pm
        
        print(f"  水{i+1}: O=({o_pos[0]:.3f}, {o_pos[1]:.3f}, {o_pos[2]:.3f})")
        print(f"       D=({d_pos[0]:.3f}, {d_pos[1]:.3f}, {d_pos[2]:.3f})")
        print(f"       位移={disp_mag:.2f} pm")
    
    print(f"\n系统能量: {energy}")
    
    # 测试能量最小化
    print("\n2. 测试能量最小化")
    print("运行能量最小化...")
    
    simulation.minimizeEnergy(maxIterations=100)
    
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    positions_minimized = state.getPositions(asNumpy=True)
    energy_min = state.getPotentialEnergy()
    
    print("\n最小化后Drude位置:")
    for i in range(2):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        o_pos = positions_minimized[o_idx].value_in_unit(unit.nanometer)
        d_pos = positions_minimized[d_idx].value_in_unit(unit.nanometer)
        
        disp = d_pos - o_pos
        disp_mag = np.linalg.norm(disp) * 1000  # pm
        
        # 计算诱导偶极矩
        q_drude = -1.71636  # e
        dipole = disp * q_drude * 4.80321  # Debye
        dipole_mag = np.linalg.norm(dipole)
        
        print(f"  水{i+1}: 位移={disp_mag:.2f} pm, 偶极矩={dipole_mag:.3f} D")
    
    print(f"\n最小化后能量: {energy_min}")
    
    # 探索其他方法
    print("\n3. 探索其他获取Drude信息的方法")
    
    # 检查DrudeForce的方法
    print("\nDrudeForce可用方法:")
    methods = [method for method in dir(drudeForce) if not method.startswith('_')]
    for method in sorted(methods):
        if 'get' in method.lower() or 'particle' in method.lower():
            print(f"  - {method}")
    
    # 获取Drude粒子信息
    print("\n获取Drude粒子参数:")
    n_particles = drudeForce.getNumParticles()
    print(f"  Drude粒子数: {n_particles}")
    
    for i in range(n_particles):
        params = drudeForce.getParticleParameters(i)
        print(f"\n  粒子{i}:")
        print(f"    Drude索引: {params[0]}")
        print(f"    Parent索引: {params[1]}")
        print(f"    电荷: {params[5]}")
        print(f"    极化率: {params[6]}")
    
    # 总结
    print("\n" + "="*70)
    print("总结")
    print("="*70)
    print("1. OpenMM可以输出Drude粒子位置")
    print("2. Drude粒子像普通粒子一样在positions数组中")
    print("3. 可以通过索引访问Drude位置并计算位移")
    print("4. DrudeLangevinIntegrator会自动更新Drude位置")
    print("5. 能量最小化也会优化Drude位置")

def test_drude_scf_mode():
    """
    测试OpenMM的SCF模式（如果有）
    """
    print("\n\n" + "="*70)
    print("探索OpenMM的SCF模式")
    print("="*70)
    
    # 检查是否有SCF相关功能
    print("\n检查OpenMM中的SCF相关功能:")
    
    # 查看DrudeForce是否有SCF相关方法
    system, topology, drudeForce = create_simple_drude_system()
    
    print("\nDrudeForce中包含'scf'的方法:")
    for attr in dir(drudeForce):
        if 'scf' in attr.lower():
            print(f"  - {attr}")
    
    # 检查是否有特殊的积分器
    print("\n检查可用的积分器:")
    integrators = []
    for name in dir(mm):
        if 'Integrator' in name and 'Drude' in name:
            integrators.append(name)
    
    for integrator in sorted(integrators):
        print(f"  - {integrator}")
    
    # 查看DrudeSCFIntegrator（如果存在）
    if hasattr(mm, 'DrudeSCFIntegrator'):
        print("\n发现DrudeSCFIntegrator!")
        scf_integrator = mm.DrudeSCFIntegrator(0.001 * unit.picosecond)
        print(f"  SCF积分器可用")
        
        # 测试SCF模式
        simulation = app.Simulation(topology, system, scf_integrator)
        positions = [[0.5, 0.5, 0.5]] * 10 * unit.nanometer
        simulation.context.setPositions(positions)
        
        print("\n使用SCF积分器优化Drude位置...")
        simulation.step(1)  # SCF只需要一步
        
        state = simulation.context.getState(getPositions=True)
        positions_scf = state.getPositions(asNumpy=True)
        
        print("SCF优化完成")
    else:
        print("\n未找到DrudeSCFIntegrator")
        print("OpenMM可能使用动力学方法而非SCF来处理Drude粒子")

def main():
    """
    主函数
    """
    print("探索OpenMM Drude粒子功能")
    print("="*70)
    
    # 测试基本的Drude位置获取
    test_drude_positions()
    
    # 探索SCF模式
    test_drude_scf_mode()

if __name__ == "__main__":
    main()