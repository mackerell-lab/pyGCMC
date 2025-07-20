#!/usr/bin/env python3
"""
分析OpenMM和PyGCMC的能量组成
理解为什么能量差异巨大但Drude位移相似
"""

import numpy as np

try:
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    print("警告：OpenMM未安装")

def analyze_openmm_energy():
    """
    分析OpenMM的能量组成
    """
    print("OpenMM能量分析")
    print("="*70)
    
    # 创建简单的2水系统
    system = mm.System()
    
    # 添加粒子
    masses = [15.99943, 0.4, 1.007947, 1.007947, 0.0]  # O, D, H1, H2, M
    for i in range(10):  # 2水×5粒子
        system.addParticle(masses[i % 5] * unit.amu)
    
    # 设置盒子
    box_size = 2.0
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0] * unit.nanometer,
        [0, box_size, 0] * unit.nanometer,
        [0, 0, box_size] * unit.nanometer
    )
    
    # 1. 只添加DrudeForce
    print("\n1. 只有DrudeForce的系统:")
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
    
    # 创建位置
    positions = [
        [0.5, 0.5, 0.5],   # O1
        [0.501, 0.5, 0.5], # D1 (稍微偏离)
        [0.596, 0.5, 0.5], # H1
        [0.452, 0.577, 0.5], # H2
        [0.5, 0.5, 0.5],   # M1
        [1.0, 1.0, 1.0],   # O2
        [1.001, 1.0, 1.0], # D2 (稍微偏离)
        [1.096, 1.0, 1.0], # H1
        [0.952, 1.077, 1.0], # H2
        [1.0, 1.0, 1.0]    # M2
    ] * unit.nanometer
    
    # 测试能量
    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(system, integrator)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy()
    forces = state.getForces()
    
    print(f"  总能量: {energy}")
    print(f"  Drude1受力: {forces[1]}")
    print(f"  Drude2受力: {forces[6]}")
    
    # 2. 添加NonbondedForce
    print("\n2. 添加NonbondedForce后:")
    
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(0.9 * unit.nanometer)
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
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
    
    # 创建新的积分器和context
    integrator2 = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(system, integrator2)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    energy_total = state.getPotentialEnergy()
    
    print(f"  总能量: {energy_total}")
    
    # 分别获取各个力的贡献
    print("\n3. 各个力的能量贡献:")
    
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force_name = force.__class__.__name__
        
        # 创建只有这个力的系统
        temp_system = mm.System()
        for j in range(10):
            temp_system.addParticle(masses[j % 5] * unit.amu)
        temp_system.setDefaultPeriodicBoxVectors(
            [box_size, 0, 0] * unit.nanometer,
            [0, box_size, 0] * unit.nanometer,
            [0, 0, box_size] * unit.nanometer
        )
        
        # 复制力（简化处理）
        if isinstance(force, mm.DrudeForce):
            temp_force = mm.DrudeForce()
            for j in range(2):
                temp_force.addParticle(j*5+1, j*5, -1, -1, -1,
                    -1.71636 * unit.elementary_charge,
                    0.0009782237 * unit.nanometer**3, 1.0, 1.0)
            temp_force.addScreenedPair(0, 1, 1.3)
        elif isinstance(force, mm.NonbondedForce):
            temp_force = mm.NonbondedForce()
            temp_force.setNonbondedMethod(mm.NonbondedForce.PME)
            temp_force.setCutoffDistance(0.9 * unit.nanometer)
            for j in range(10):
                charge = charges[j % 5] * unit.elementary_charge
                sigma = 0.318395 * unit.nanometer if j % 5 == 0 else 1.0 * unit.nanometer
                epsilon = 0.88257 * unit.kilojoule_per_mole if j % 5 == 0 else 0.0 * unit.kilojoule_per_mole
                temp_force.addParticle(charge, sigma, epsilon)
            for j in range(2):
                base = j * 5
                for k in range(5):
                    for l in range(k+1, 5):
                        temp_force.addException(base+k, base+l, 0, 1, 0)
        else:
            continue
        
        temp_system.addForce(temp_force)
        
        temp_integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
        temp_context = mm.Context(temp_system, temp_integrator)
        temp_context.setPositions(positions)
        
        temp_state = temp_context.getState(getEnergy=True)
        temp_energy = temp_state.getPotentialEnergy()
        
        print(f"  {force_name}: {temp_energy}")
    
    # 使用SCF积分器
    print("\n4. 使用DrudeSCFIntegrator:")
    
    scf_integrator = mm.DrudeSCFIntegrator(0.001 * unit.picosecond)
    scf_context = mm.Context(system, scf_integrator)
    
    # 重置Drude到parent位置
    positions_reset = list(positions)
    positions_reset[1] = positions[0]  # D1 = O1
    positions_reset[6] = positions[5]  # D2 = O2
    
    scf_context.setPositions(positions_reset)
    
    # 初始能量
    state_init = scf_context.getState(getEnergy=True)
    energy_init = state_init.getPotentialEnergy()
    print(f"  初始能量（Drude在parent）: {energy_init}")
    
    # 运行SCF
    scf_integrator.step(1)
    
    state_final = scf_context.getState(getPositions=True, getEnergy=True)
    energy_final = state_final.getPotentialEnergy()
    positions_final = state_final.getPositions()
    
    print(f"  SCF后能量: {energy_final}")
    print(f"  能量变化: {energy_final - energy_init}")
    
    # 检查Drude位移
    pos1_o = positions_final[0].value_in_unit(unit.nanometer)
    pos1_d = positions_final[1].value_in_unit(unit.nanometer)
    pos2_o = positions_final[5].value_in_unit(unit.nanometer)
    pos2_d = positions_final[6].value_in_unit(unit.nanometer)
    
    d1_disp = np.linalg.norm(pos1_d - pos1_o) * 10  # Angstrom
    d2_disp = np.linalg.norm(pos2_d - pos2_o) * 10
    
    print(f"  Drude1位移: {d1_disp:.3f} Å")
    print(f"  Drude2位移: {d2_disp:.3f} Å")

def analyze_energy_difference():
    """
    分析能量差异的原因
    """
    print("\n\n能量差异分析")
    print("="*70)
    
    print("\n可能的原因：")
    print("1. OpenMM的DrudeSCFIntegrator可能只优化Drude位置，不计算完整能量")
    print("2. PyGCMC包含了所有能量项（LJ + 库仑 + Drude）")
    print("3. 单位或常数差异")
    print("4. 截断处理不同")
    
    print("\n验证方法：")
    print("1. 分别计算各个能量项")
    print("2. 检查是否包含了所有相互作用")
    print("3. 比较相同配置下的能量")

def main():
    """
    主函数
    """
    print("OpenMM和PyGCMC能量差异分析")
    print("="*70)
    
    if HAS_OPENMM:
        analyze_openmm_energy()
    
    analyze_energy_difference()

if __name__ == "__main__":
    main()