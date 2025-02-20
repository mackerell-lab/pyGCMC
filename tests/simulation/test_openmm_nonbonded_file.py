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
    
    # 打印系统基本信息（只打印一次）
    print(f"\n系统信息: {system.getNumParticles()}个粒子")
    
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

        # 2. 计算自定义公式的能量
        # 2.1 计算正常非键相互作用
        combined_nonbonded_expression = """
            (kC * q1 * q2 / r + 4 * sqrt(eps1*eps2) * ((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6))
        """

        nonbonded_force = CustomNonbondedForce(combined_nonbonded_expression.replace('\n', '').strip())
        nonbonded_force.addPerParticleParameter("q")
        nonbonded_force.addPerParticleParameter("sigma")
        nonbonded_force.addPerParticleParameter("eps")
        nonbonded_force.addGlobalParameter("kC", 138.935456)
        nonbonded_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)

        # 添加粒子参数
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nonbonded_force.addParticle([charge, sigma, epsilon])

        # 2.2 计算例外相互作用
        exception_expression = """
            (kC * chargeprod / r * coulombscale + 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * ljscale)
        """
        
        exception_force = CustomBondForce(exception_expression)
        exception_force.addPerBondParameter("chargeprod")
        exception_force.addPerBondParameter("sigma")
        exception_force.addPerBondParameter("epsilon")
        exception_force.addPerBondParameter("ljscale")
        exception_force.addPerBondParameter("coulombscale")
        exception_force.addGlobalParameter("kC", 138.935456)

        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            coulomb_scale = 1.0
            lj_scale = 1.0
            exception_force.addBond(p1, p2, [chargeProd, sigma, epsilon, lj_scale, coulomb_scale])
            nonbonded_force.addExclusion(p1, p2)

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
        
        print(f"{dist:6.2f}  {openmm_val:10.4f}  {custom_val:11.4f}  {rel_error:8.4f}")
        
        # 验证结果
        assert rel_error < 1e-4, f"Relative error too large: {rel_error:.4f}%"

