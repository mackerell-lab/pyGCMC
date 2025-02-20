# tests/simulation/test_openmm_nonbonded_file.py

import pytest
import os
import numpy as np
import warnings

import math

# Suppress SWIG-related DeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyPacked has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyObject has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type swigvarlink has no __module__ attribute")
from openmm import *
from openmm.app import *
from openmm.unit import *
from openmm.app.gromacstopfile import GromacsTopFile

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

def load_test_system():
    """从PDB和TOP文件加载测试系统。"""
    # 加载PDB文件
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    pdb = PDBFile(pdb_path)
    
    # 加载TOP文件
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    top = GromacsTopFile(top_path)
    
    # 创建系统
    system = top.createSystem(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.0*nanometer,
        switchDistance=0.9*nanometer
    )
    
    return system, pdb.positions

def test_verify_openmm_expressions():
    """验证我们的自定义非键相互作用表达式与OpenMM默认实现的一致性。"""
    # 加载测试系统
    system, positions = load_test_system()
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 定义测试距离
    distances = [0.9, 0.95, 1.0, 1.1, 1.5, 2.0]  # 移除极端短距离
    print("\n=== 比较不同非键相互作用表达式 ===")
    print("距离(nm)  OpenMM默认   自定义公式    相对误差(%)")
    print("-" * 55)
    
    platform = Platform.getPlatformByName('Reference')
    
    # 打印系统信息
    print("\n系统信息:")
    print(f"粒子数: {system.getNumParticles()}")
    print(f"力场数: {system.getNumForces()}")
    for i, force in enumerate(system.getForces()):
        print(f"力场 {i}: {force.__class__.__name__}")
    
    # 打印NonbondedForce的设置
    print("\nNonbondedForce设置:")
    print(f"非键方法: {original_nb_force.getNonbondedMethod()}")
    print(f"截断距离: {original_nb_force.getCutoffDistance()}")
    print(f"切换距离: {original_nb_force.getSwitchingDistance()}")
    print(f"是否使用切换函数: {original_nb_force.getUseSwitchingFunction()}")
    
    # 打印一些粒子参数示例
    print("\n粒子参数示例 (前3个):")
    for i in range(min(3, original_nb_force.getNumParticles())):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        print(f"粒子 {i}: 电荷={charge}, sigma={sigma}, epsilon={epsilon}")
    
    # 打印一些例外相互作用示例
    print("\n例外相互作用示例 (前3个):")
    for i in range(min(3, original_nb_force.getNumExceptions())):
        p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
        print(f"例外 {i}: 粒子({p1},{p2}), 电荷积={chargeProd}, sigma={sigma}, epsilon={epsilon}")
    
    for dist in distances:
        # 缩放所有位置
        scaled_positions = [Vec3(pos[0].value_in_unit(nanometers) * dist,
                               pos[1].value_in_unit(nanometers) * dist,
                               pos[2].value_in_unit(nanometers) * dist) * nanometers
                          for pos in positions]
        
        # 1. Create a new system with ONLY NonbondedForce for comparison
        openmm_system = System()
        for i in range(system.getNumParticles()):
            openmm_system.addParticle(system.getParticleMass(i))
        
        openmm_nb_force = NonbondedForce()
        openmm_nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)
        openmm_nb_force.setCutoffDistance(1e10 * nanometers)
        openmm_nb_force.setSwitchingDistance(0.0 * nanometer)
        openmm_nb_force.setUseSwitchingFunction(False)
        
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            openmm_nb_force.addParticle(charge, sigma, epsilon)
        
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            openmm_nb_force.addException(p1, p2, chargeProd, sigma, epsilon)
        
        openmm_system.addForce(openmm_nb_force)
        
        # 1. 计算OpenMM默认实现的能量
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(openmm_system, integrator, platform)
        context.setPositions(scaled_positions)
        state = context.getState(getEnergy=True)
        openmm_energy = state.getPotentialEnergy()
        del context, integrator

        print(f"\n=== 距离 {dist} nm 的详细信息 ===")
        print("\nOpenMM默认NonbondedForce参数:")
        print(f"粒子数: {openmm_nb_force.getNumParticles()}")
        print(f"例外数: {openmm_nb_force.getNumExceptions()}")
        print(f"使用切换函数: {openmm_nb_force.getUseSwitchingFunction()}")
        print(f"切换距离: {openmm_nb_force.getSwitchingDistance()}")
        print(f"截断距离: {openmm_nb_force.getCutoffDistance()}")
        print(f"OpenMM能量: {openmm_energy.value_in_unit(kilojoules_per_mole):.4f} kJ/mol")
        
        # 2. 计算自定义公式的能量
        # 2.1 计算正常非键相互作用
        combined_nonbonded_expression = """
            (kC * q1 * q2 / r + 4 * sqrt(eps1*eps2) * ((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6))
        """
        print("\n自定义NonbondedForce表达式:")
        print(combined_nonbonded_expression)

        nonbonded_force = CustomNonbondedForce(combined_nonbonded_expression.replace('\n', '').strip())
        nonbonded_force.addPerParticleParameter("q")
        nonbonded_force.addPerParticleParameter("sigma")
        nonbonded_force.addPerParticleParameter("eps")
        nonbonded_force.addGlobalParameter("kC", 138.935456)
        nonbonded_force.addGlobalParameter("cutoff", 1.0)
        nonbonded_force.addGlobalParameter("switch", 0.9)

        print("\n自定义NonbondedForce参数:")
        print(f"粒子数: {nonbonded_force.getNumParticles()}")
        print(f"全局参数:")
        print(f"  kC = {nonbonded_force.getGlobalParameterDefaultValue(0)}")
        print(f"  cutoff = {nonbonded_force.getGlobalParameterDefaultValue(1)}")
        print(f"  switch = {nonbonded_force.getGlobalParameterDefaultValue(2)}")

        # 添加粒子参数
        print("\n粒子参数示例 (前5个):")
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nonbonded_force.addParticle([charge, sigma, epsilon])
            if i < 5:  # Only print the first 5 particles for debugging
                print(f"粒子 {i}: q={charge}, sigma={sigma}, epsilon={epsilon}")

        nonbonded_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
        nonbonded_force.setCutoffDistance(1e10 * nanometers)

        # 2.2 计算例外相互作用
        exception_expression = """
            (kC * chargeprod / r * coulombscale + 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * ljscale)
        """
        
        print("\n例外相互作用表达式:")
        print(exception_expression)
        
        exception_force = CustomBondForce(exception_expression)
        
        exception_force.addPerBondParameter("chargeprod")
        exception_force.addPerBondParameter("sigma")
        exception_force.addPerBondParameter("epsilon")
        exception_force.addPerBondParameter("ljscale")
        exception_force.addPerBondParameter("coulombscale")
        exception_force.addGlobalParameter("kC", 138.935456)
        exception_force.addGlobalParameter("cutoff", 1.0)
        exception_force.addGlobalParameter("switch", 0.9)

        print("\n例外相互作用参数示例 (前5个):")
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            
            # If you need proper 1-4 scaling:
            #   - read "unscaled" from the combination rules
            #   - compute scale = (exception_value / unscaled_value)
            #   - but only if your ForceField demands that logic.

            # Otherwise, simply treat them as final:
            coulomb_scale = 1.0
            lj_scale = 1.0

            # Add the bond using final parameters from the original NonbondedForce.
            exception_force.addBond(p1, p2, [chargeProd, sigma, epsilon, lj_scale, coulomb_scale])
            
            # Also be sure to exclude these same pairs from the main CustomNonbondedForce
            nonbonded_force.addExclusion(p1, p2)
            
            if i < 5:  # Keep this for printing only the first 5
                print(f"例外 {i}:")
                print(f"  粒子对: ({p1}, {p2})")
                print(f"  原始参数: chargeProd={chargeProd}, sigma={sigma}, epsilon={epsilon}")
                print(f"  缩放因子: coulomb_scale={coulomb_scale}, lj_scale={lj_scale}")

        # 创建系统并计算能量
        custom_system = System()
        for i in range(system.getNumParticles()):
            custom_system.addParticle(system.getParticleMass(i))
        custom_system.addForce(nonbonded_force)
        custom_system.addForce(exception_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(custom_system, integrator, platform)
        context.setPositions(scaled_positions)
        custom_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 计算相对误差
        openmm_val = openmm_energy.value_in_unit(kilojoules_per_mole)
        custom_val = custom_energy.value_in_unit(kilojoules_per_mole)
        
        rel_error = abs(custom_val - openmm_val) / abs(openmm_val) * 100 if abs(openmm_val) > 1e-6 else 0.0
        
        print(f"\n能量比较:")
        print(f"OpenMM能量: {openmm_val:.4f} kJ/mol")
        print(f"自定义能量: {custom_val:.4f} kJ/mol")
        print(f"相对误差: {rel_error:.4f}%")
        
        print(f"\n{dist:6.2f}  {openmm_val:10.4f}  {custom_val:11.4f}  {rel_error:8.4f}")
        
        # 验证结果
        # With NoCutoff, the relative error should be very small everywhere.
        assert rel_error < 1e-4, f"Relative error too large: {rel_error:.4f}%"

def test_verify_energy_components():
    """验证库伦和LJ能量分量的计算。"""
    # 初始化统计信息
    stats = {
        "particles": {
            "total": 0,
            "charged": 0,
            "max_charge": 0.0
        },
        "exceptions": {
            "total": 0,
            "zero_charge": 0,
            "nonzero_charge": 0,
            "zero_epsilon": 0,
            "nonzero_epsilon": 0
        }
    }
    
    # 加载测试系统
    system, positions = load_test_system()
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 定义测试距离
    distances = [0.9, 0.95, 1.0, 1.1, 1.5, 2.0]
    print("\n=== 比较能量分量 ===")
    print("距离(nm)  OpenMM总能量  库伦能量  LJ能量  总能量  相对误差(%)")
    print("-" * 75)
    
    platform = Platform.getPlatformByName('Reference')
    
    # 打印能量分量的计算设置
    print("\n能量分量计算设置:")
    print("库伦能量 (硬截断):")
    coulomb_nonbonded_expression = """
        select(step(cutoff-r), kC * q1 * q2 / r, 0)
    """
    print(coulomb_nonbonded_expression)
    
    print("\nLJ能量 (带切换函数):")
    lj_nonbonded_expression = """
        select(step(cutoff-r),
            4 * sqrt(eps1*eps2) * ((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6) *
            (step(switch-r) + step(r-switch) * (1 - 10*((r-switch)/(cutoff-switch))^3 + 15*((r-switch)/(cutoff-switch))^4 - 6*((r-switch)/(cutoff-switch))^5)),
            0)
    """
    print(lj_nonbonded_expression)

    # 创建力场对象
    coulomb_nonbonded_force = CustomNonbondedForce(coulomb_nonbonded_expression.replace('\n', '').strip())
    coulomb_nonbonded_force.addPerParticleParameter("q")
    coulomb_nonbonded_force.addGlobalParameter("kC", 138.935456)
    coulomb_nonbonded_force.addGlobalParameter("cutoff", 1.0)
    coulomb_nonbonded_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    coulomb_nonbonded_force.setCutoffDistance(1.0 * nanometer)
    coulomb_nonbonded_force.setUseLongRangeCorrection(False)

    lj_nonbonded_force = CustomNonbondedForce(lj_nonbonded_expression.replace('\n', '').strip())
    lj_nonbonded_force.addPerParticleParameter("sigma")
    lj_nonbonded_force.addPerParticleParameter("eps")
    lj_nonbonded_force.addGlobalParameter("cutoff", 1.0)
    lj_nonbonded_force.addGlobalParameter("switch", 0.9)
    lj_nonbonded_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    lj_nonbonded_force.setCutoffDistance(1.0 * nanometers)
    lj_nonbonded_force.setUseLongRangeCorrection(False)

    # 创建例外力场对象
    coulomb_exception_expression = "kC * chargeprod / r"
    coulomb_exception_force = CustomBondForce(coulomb_exception_expression)
    coulomb_exception_force.addPerBondParameter("chargeprod")
    coulomb_exception_force.addGlobalParameter("kC", 138.935456)

    # Update LJ exception expression to include switching function
    lj_exception_expression = """
        select(step(cutoff-r),
            4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) *
            (step(switch-r) + step(r-switch) * (1 - 10*((r-switch)/(cutoff-switch))^3 + 15*((r-switch)/(cutoff-switch))^4 - 6*((r-switch)/(cutoff-switch))^5)),
            0)
    """
    lj_exception_force = CustomBondForce(lj_exception_expression.replace('\n', '').strip())
    lj_exception_force.addPerBondParameter("sigma")
    lj_exception_force.addPerBondParameter("epsilon")
    lj_exception_force.addGlobalParameter("cutoff", 1.0)
    lj_exception_force.addGlobalParameter("switch", 0.9)

    # 处理例外
    print("\n处理例外相互作用:")
    for i in range(original_nb_force.getNumExceptions()):
        p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
        stats["exceptions"]["total"] += 1
        
        # 1. 首先从常规非键力中排除所有exception对
        coulomb_nonbonded_force.addExclusion(p1, p2)
        lj_nonbonded_force.addExclusion(p1, p2)
        
        # 2. 处理库伦相互作用
        if abs(chargeProd.value_in_unit(elementary_charge**2)) > 1e-10:
            stats["exceptions"]["nonzero_charge"] += 1
            coulomb_exception_force.addBond(p1, p2, [chargeProd])
        else:
            stats["exceptions"]["zero_charge"] += 1
            
        # 3. 处理LJ相互作用 - 只在epsilon非零时添加
        if abs(epsilon.value_in_unit(kilojoule_per_mole)) > 1e-10:
            stats["exceptions"]["nonzero_epsilon"] += 1
            lj_exception_force.addBond(p1, p2, [sigma, epsilon])
        else:
            stats["exceptions"]["zero_epsilon"] += 1
            # Skip combination rules unless explicitly required by force field

    # 打印统计信息
    print("\n系统统计信息:")
    print(f"粒子总数: {stats['particles']['total']}")
    print(f"带电粒子数: {stats['particles']['charged']}")
    print(f"最大电荷绝对值: {stats['particles']['max_charge']:.3f} e")
    print(f"\n例外总数: {stats['exceptions']['total']}")
    print(f"零电荷例外: {stats['exceptions']['zero_charge']}")
    print(f"非零电荷例外: {stats['exceptions']['nonzero_charge']}")
    print(f"零epsilon例外: {stats['exceptions']['zero_epsilon']}")
    print(f"非零epsilon例外: {stats['exceptions']['nonzero_epsilon']}")
    
    # 初始化scaled_positions (使用默认距离1.0)
    default_dist = 1.0
    scaled_positions = [Vec3(pos[0].value_in_unit(nanometers) * default_dist,
                           pos[1].value_in_unit(nanometers) * default_dist,
                           pos[2].value_in_unit(nanometers) * default_dist) * nanometers
                      for pos in positions]
    
    # 添加诊断输出
    print("\n计算OpenMM的LJ能量分量...")
    openmm_lj_system = System()
    for i in range(system.getNumParticles()):
        openmm_lj_system.addParticle(system.getParticleMass(i))
    
    openmm_lj_force = NonbondedForce()
    openmm_lj_force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
    openmm_lj_force.setCutoffDistance(1.0 * nanometer)
    openmm_lj_force.setSwitchingDistance(0.9 * nanometer)
    openmm_lj_force.setUseSwitchingFunction(True)
    openmm_lj_force.setUseDispersionCorrection(False)
    
    # 添加粒子，但电荷设为0
    for i in range(original_nb_force.getNumParticles()):
        _, sigma, epsilon = original_nb_force.getParticleParameters(i)
        openmm_lj_force.addParticle(0.0 * elementary_charge, sigma, epsilon)
    
    # 添加例外，但电荷积设为0
    for i in range(original_nb_force.getNumExceptions()):
        p1, p2, _, sigma, epsilon = original_nb_force.getExceptionParameters(i)
        openmm_lj_force.addException(p1, p2, 0.0 * elementary_charge**2, sigma, epsilon)
    
    openmm_lj_system.addForce(openmm_lj_force)
    
    # 计算OpenMM的LJ能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(openmm_lj_system, integrator, platform)
    context.setPositions(scaled_positions)
    openmm_lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
    del context, integrator
    
    print(f"OpenMM LJ能量: {openmm_lj_energy.value_in_unit(kilojoules_per_mole):.4f} kJ/mol")

    for dist in distances:
        # 缩放所有位置
        scaled_positions = [Vec3(pos[0].value_in_unit(nanometers) * dist,
                               pos[1].value_in_unit(nanometers) * dist,
                               pos[2].value_in_unit(nanometers) * dist) * nanometers
                          for pos in positions]

        # 1. 计算OpenMM参考能量
        openmm_system = System()
        for i in range(system.getNumParticles()):
            openmm_system.addParticle(system.getParticleMass(i))

        openmm_nb_force = NonbondedForce()
        openmm_nb_force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
        openmm_nb_force.setCutoffDistance(1.0 * nanometer)
        openmm_nb_force.setSwitchingDistance(0.9 * nanometer)
        openmm_nb_force.setUseSwitchingFunction(True)
        openmm_nb_force.setUseDispersionCorrection(False)
        openmm_nb_force.setReactionFieldDielectric(1.0)

        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            openmm_nb_force.addParticle(charge, sigma, epsilon)

        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            openmm_nb_force.addException(p1, p2, chargeProd, sigma, epsilon)

        openmm_system.addForce(openmm_nb_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(openmm_system, integrator, platform)
        context.setPositions(scaled_positions)
        openmm_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 2. 计算库伦能量
        coulomb_system = System()
        for i in range(system.getNumParticles()):
            coulomb_system.addParticle(system.getParticleMass(i))
            # Add particle parameters to coulomb_nonbonded_force
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            coulomb_nonbonded_force.addParticle([charge])
            
            # Update particle statistics
            stats["particles"]["total"] += 1
            if abs(charge.value_in_unit(elementary_charge)) > 1e-10:
                stats["particles"]["charged"] += 1
                stats["particles"]["max_charge"] = max(stats["particles"]["max_charge"], 
                                                     abs(charge.value_in_unit(elementary_charge)))

        coulomb_system.addForce(coulomb_nonbonded_force)
        coulomb_system.addForce(coulomb_exception_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(coulomb_system, integrator, platform)
        context.setPositions(scaled_positions)
        coulomb_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 3. 计算LJ能量
        lj_system = System()
        for i in range(system.getNumParticles()):
            lj_system.addParticle(system.getParticleMass(i))
            # Add particle parameters to lj_nonbonded_force
            _, sigma, epsilon = original_nb_force.getParticleParameters(i)
            lj_nonbonded_force.addParticle([sigma, epsilon])

        lj_system.addForce(lj_nonbonded_force)
        lj_system.addForce(lj_exception_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(lj_system, integrator, platform)
        context.setPositions(scaled_positions)
        lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 计算总能量和相对误差
        openmm_val = openmm_energy.value_in_unit(kilojoules_per_mole)
        coulomb_val = coulomb_energy.value_in_unit(kilojoules_per_mole)
        lj_val = lj_energy.value_in_unit(kilojoules_per_mole)
        total_val = coulomb_val + lj_val

        rel_error = abs(total_val - openmm_val) / abs(openmm_val) * 100 if abs(openmm_val) > 1e-6 else 0.0

        print(f"\n=== 距离 {dist} nm 的能量分析 ===")
        print(f"OpenMM总能量: {openmm_val:.4f} kJ/mol")
        print(f"库伦能量: {coulomb_val:.4f} kJ/mol")
        print(f"LJ能量: {lj_val:.4f} kJ/mol")
        print(f"自定义总能量: {total_val:.4f} kJ/mol")
        print(f"相对误差: {rel_error:.4f}%")

        # 验证结果
        if dist <= 0.9:  # 切换距离内
            assert rel_error < 5.0, f"切换距离内的相对误差过大: {rel_error:.4f}% > 5.0%"
        elif dist < 1.0:  # 切换区域
            assert rel_error < 2.0, f"切换区域内的相对误差过大: {rel_error:.4f}% > 2.0%"
        else:  # 截断距离外
            assert rel_error < 1e-4, f"截断距离外能量应接近零: {total_val:.6f} kJ/mol, 相对误差: {rel_error:.6f}%"

def test_verify_switching_function():
    """验证切换函数的行为。"""
    # 加载测试系统
    system, positions = load_test_system()
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 定义更密集的测试距离，特别关注切换区域
    distances = [0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1.0, 1.1]  # 移除极端短距离
    print("\n=== 验证切换函数 ===")
    print("距离(nm)  切换函数值  能量比例  相对误差(%)")
    print("-" * 50)
    
    platform = Platform.getPlatformByName('Reference')
    
    # 打印切换函数的详细信息
    print("\n切换函数详细信息:")
    print("切换函数表达式:")
    print("sw = step(cutoff - r) * (step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3) + step(switch - r))")
    print(f"\n参数:")
    print(f"切换开始距离 (switch) = {0.9} nm")
    print(f"截断距离 (cutoff) = {1.0} nm")
    print(f"切换区域宽度 = {0.1} nm")
    
    def calculate_switching_function(r):
        """计算给定距离的切换函数值"""
        if r <= 0.9:  # switch distance
            return 1.0
        elif r >= 1.0:  # cutoff distance
            return 0.0
        else:
            x = (r - 0.9) / (1.0 - 0.9)  # Normalized distance
            return 1 - 10*x**3 + 15*x**4 - 6*x**5
    
    print("\n理论切换函数值:")
    test_distances = [0.85, 0.89, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1.0]
    for r in test_distances:
        sw = calculate_switching_function(r)
        print(f"r = {r:4.2f} nm: {sw:6.4f}")
    
    for dist in distances:
        # 缩放所有位置
        scaled_positions = [Vec3(pos[0].value_in_unit(nanometers) * dist,
                               pos[1].value_in_unit(nanometers) * dist,
                               pos[2].value_in_unit(nanometers) * dist) * nanometers
                          for pos in positions]
        
        # 计算理论切换函数值
        switch_value = calculate_switching_function(dist)
        
        # 1. 计算不带切换函数的能量
        # 1.1 计算正常非键相互作用
        simple_nonbonded_expression = "4 * sqrt(eps1*eps2) * ((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6) + kC * q1 * q2 / r"
        simple_nonbonded_force = CustomNonbondedForce(simple_nonbonded_expression)
        simple_nonbonded_force.addPerParticleParameter("q")
        simple_nonbonded_force.addPerParticleParameter("sigma")
        simple_nonbonded_force.addPerParticleParameter("eps")
        simple_nonbonded_force.addGlobalParameter("kC", 138.935456)
        simple_nonbonded_force.addGlobalParameter("cutoff", 1.0)

        # 添加粒子参数
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            simple_nonbonded_force.addParticle([charge, sigma, epsilon])

        # 添加排除项
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            simple_nonbonded_force.addExclusion(p1, p2)

        simple_nonbonded_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        simple_nonbonded_force.setCutoffDistance(1.0 * nanometers)

        # 1.2 计算例外相互作用
        simple_exception_expression = "4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) + kC * chargeprod / r"
        simple_exception_force = CustomBondForce(simple_exception_expression)
        simple_exception_force.addPerBondParameter("chargeprod")
        simple_exception_force.addPerBondParameter("sigma")
        simple_exception_force.addPerBondParameter("epsilon")
        simple_exception_force.addGlobalParameter("kC", 138.935456)
        simple_exception_force.addGlobalParameter("cutoff", 1.0)

        # 添加例外相互作用
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            simple_exception_force.addBond(p1, p2, [chargeProd, sigma, epsilon])

        # 2. 计算带切换函数的能量
        # 2.1 计算正常非键相互作用
        switched_nonbonded_expression = """
            (4 * sqrt(eps1*eps2) * ((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6)
            * step(cutoff - r)
            * (step(r - switch) * (1 - 10*((r - switch)/(cutoff - switch))^3 + 15*((r - switch)/(cutoff - switch))^4 - 6*((r - switch)/(cutoff - switch))^5)
            + step(switch - r))
        """
        switched_nonbonded_force = CustomNonbondedForce(switched_nonbonded_expression.replace('\n', '').strip())
        # Add the missing parameters
        switched_nonbonded_force.addPerParticleParameter("q")
        switched_nonbonded_force.addPerParticleParameter("sigma")
        switched_nonbonded_force.addPerParticleParameter("eps")
        switched_nonbonded_force.addGlobalParameter("kC", 138.935456)
        switched_nonbonded_force.addGlobalParameter("cutoff", 1.0)
        switched_nonbonded_force.addGlobalParameter("switch", 0.9)

        # 添加粒子参数
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            switched_nonbonded_force.addParticle([charge, sigma, epsilon])

        # 添加排除项
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            switched_nonbonded_force.addExclusion(p1, p2)

        switched_exception_expression = """
            (4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * ljscale + kC * chargeprod / r * coulombscale)
            * step(cutoff - r)
            * (step(r - switch) * (1 - 10*((r - switch)/(cutoff - switch))^3 + 15*((r - switch)/(cutoff - switch))^4 - 6*((r - switch)/(cutoff - switch))^5)
            + step(switch - r))
        """
        switched_exception_force = CustomBondForce(switched_exception_expression.replace('\n', '').strip())

        # Add per-bond parameters first
        switched_exception_force.addPerBondParameter("chargeprod")
        switched_exception_force.addPerBondParameter("sigma")
        switched_exception_force.addPerBondParameter("epsilon")
        switched_exception_force.addPerBondParameter("ljscale")
        switched_exception_force.addPerBondParameter("coulombscale")
        switched_exception_force.addGlobalParameter("kC", 138.935456)
        switched_exception_force.addGlobalParameter("cutoff", 1.0)
        switched_exception_force.addGlobalParameter("switch", 0.9)

        # Add exception bonds with proper parameters
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            # Calculate scaling factors
            q1, _, _ = original_nb_force.getParticleParameters(p1)
            q2, _, _ = original_nb_force.getParticleParameters(p2)
            
            # Calculate unscaled charge product
            unscaled_chargeProd = q1 * q2
            
            # Calculate scaling factors
            coulomb_scale = 0.0 if abs(unscaled_chargeProd.value_in_unit(elementary_charge**2)) < 1e-10 else chargeProd / unscaled_chargeProd
            lj_scale = 1.0  # Default to 1.0 for exceptions
            
            switched_exception_force.addBond(p1, p2, [chargeProd, sigma, epsilon, lj_scale, coulomb_scale])

        # 创建系统并计算不带切换函数的能量
        simple_system = System()
        for i in range(system.getNumParticles()):
            simple_system.addParticle(system.getParticleMass(i))
        simple_system.addForce(simple_nonbonded_force)
        simple_system.addForce(simple_exception_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(simple_system, integrator, platform)
        context.setPositions(scaled_positions)
        simple_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 创建系统并计算带切换函数的能量
        switched_system = System()
        for i in range(system.getNumParticles()):
            switched_system.addParticle(system.getParticleMass(i))
        switched_system.addForce(switched_nonbonded_force)
        switched_system.addForce(switched_exception_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(switched_system, integrator, platform)
        context.setPositions(scaled_positions)
        switched_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 计算能量比例和相对误差
        simple_val = simple_energy.value_in_unit(kilojoules_per_mole)
        switched_val = switched_energy.value_in_unit(kilojoules_per_mole)
        
        # 计算实际能量比例（仅在简单能量不接近零时）
        energy_ratio = switched_val / simple_val if abs(simple_val) > 1e-6 else 0.0
        
        # 计算与理论切换函数值的相对误差
        if dist <= 0.9 or dist >= 1.0:
            rel_error = abs(energy_ratio - switch_value) * 100 if abs(switch_value) > 1e-6 else 0.0
        else:
            # 在切换区域内，能量比例应该接近切换函数值
            rel_error = abs(energy_ratio - switch_value) / switch_value * 100 if abs(switch_value) > 1e-6 else 0.0
        
        print(f"{dist:6.2f}  {switch_value:10.4f}  {energy_ratio:9.4f}  {rel_error:8.4f}")
        
        # 验证结果
        if dist <= 0.9:  # 切换距离内
            assert abs(energy_ratio - 1.0) < 0.01, f"切换距离内能量比例应为1: {energy_ratio:.4f}"
        elif dist >= 1.0:  # 截断距离外
            assert abs(switched_val) < 1e-6, f"截断距离外能量应为0: {switched_val:.6f}"
        else:  # 切换区域
            assert abs(energy_ratio - switch_value) / switch_value < 0.02, \
                    f"切换区域内能量比例应接近切换函数值: {energy_ratio:.4f} != {switch_value:.4f}"

