# tests/simulation/test_switching_function.py

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCInfo, MCMovementResidueInfo
import os
import math
import sys

# 设置日志级别为INFO，以便查看调试输出
pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)

def create_test_system(box_size=4.0):
    """
    Create a simple test system with two atoms
    
    Args:
        box_size: box size (nm)
    Returns:
        state: MC state object with two atoms
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff (12 Å) - 使用 CHARMM 的 ctofnb 值作为截断
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Type 0 and Type 1
    
    # LJ parameters
    sigma = 0.3  # nm
    eps = 0.5    # kJ/mol
    
    # Set LJ parameter matrix (2x2 matrix flattened to array)
    ff.ljSigma = [
        sigma, sigma,  # type 0 to types 0, 1
        sigma, sigma   # type 1 to types 0, 1
    ]
    
    ff.ljEps = [
        eps, eps,  # type 0 to types 0, 1
        eps, eps   # type 1 to types 0, 1
    ]
    
    state.forcefield = ff
    
    # First atom
    atom1 = MCAtom()
    atom1.type = 0
    atom1.charge = 0.0  # Neutral
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    
    # Second atom
    atom2 = MCAtom()
    atom2.type = 1
    atom2.charge = 0.0  # Neutral
    atom2.x = 1.0  # Will be adjusted in tests
    atom2.y = 0.0
    atom2.z = 0.0
    
    # Add atoms to state
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # Create two residues and atoms
    res1 = MCResidue()
    res1.active = True
    res1.atomStart = 0
    res1.atomCount = 1
    res1.type = 0
    
    res2 = MCResidue()
    res2.active = True
    res2.atomStart = 1
    res2.atomCount = 1
    res2.type = 1
    
    # Add residues to state
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    return state

def test_switching_function():
    """
    测试CHARMM平滑函数的计算
    """
    # 在不同的距离下计算平滑函数的值
    # 按照 CHARMM 的约定:
    # ctonnb = 1.0 - inner cutoff (r_on)，平滑衰减开始处 (10 Å)
    # ctofnb = 1.2 - outer cutoff (r_off)，势能最终衰减至 0 处 (12 Å)
    r_on = 1.0  # CHARMM 的 ctonnb (10 Å)
    r_off = 1.2  # CHARMM 的 ctofnb (12 Å)
    
    # 创建测试系统
    state = create_test_system()
    
    # 启用平滑函数
    pygcmc.enableSwitchingFunction(state, True, r_on, r_off)
    
    # 计算一系列距离上的平滑函数值
    distances = [0.9 + i * 0.4/50 for i in range(50)]  # 从0.9到1.3的50个点
    switch_values = [pygcmc.calculateSwitchingFunction(state, r) for r in distances]
    
    # 验证关键点的值
    for r, s in zip(distances, switch_values):
        if r <= r_on:
            assert s == pytest.approx(1.0), f"S({r}) should be 1.0 when r <= r_on"
        elif r >= r_off:
            assert s == pytest.approx(0.0), f"S({r}) should be 0.0 when r >= r_off"
        else:
            # 根据公式计算预期值进行验证
            r2 = r * r
            ron2 = r_on * r_on
            roff2 = r_off * r_off
            
            numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0*r2 - 3.0*ron2)
            denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2)
            expected = numerator / denominator
            
            # 增加容差值以适应浮点精度差异
            assert s == pytest.approx(expected, abs=1e-6), f"S({r}) calculation error"
            
    # 确保平滑函数在区间边界处连续（在r_on处值为1.0，在r_off处值为0.0）
    s_at_ron = pygcmc.calculateSwitchingFunction(state, r_on)
    s_at_roff = pygcmc.calculateSwitchingFunction(state, r_off)
    
    assert s_at_ron == pytest.approx(1.0), f"S({r_on}) should be 1.0"
    assert s_at_roff == pytest.approx(0.0), f"S({r_off}) should be 0.0"

def test_energy_with_switching():
    """
    测试使用平滑截断计算能量
    """
    # 创建测试系统
    state = create_test_system()
    
    # 设置平滑函数参数
    # 按照 CHARMM 的约定:
    # ctonnb = 1.0 - inner cutoff (r_on)，平滑衰减开始处 (10 Å)
    # ctofnb = 1.2 - outer cutoff (r_off)，势能最终衰减至 0 处 (12 Å)
    r_on = 1.0   # 内截断半径 (10 Å)
    r_off = 1.2  # 外截断半径 (12 Å)
    
    # 设置第二个原子的位置
    distances = [0.9 + i * 0.4/50 for i in range(50)]  # 从0.9到1.3的50个点
    
    # 计算标准硬截断LJ能量
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)  # 禁用平滑函数
    
    energies_hard_cutoff = []
    for r in distances:
        state.atoms[1].x = r  # 设置距离
        pygcmc.computeSystemEnergy(state)  # 使用默认方法计算能量
        energies_hard_cutoff.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 计算带平滑函数的LJ能量
    pygcmc.enableSwitchingFunction(state, True, r_on, r_off)  # 启用平滑函数
    
    energies_with_switching = []
    for r in distances:
        state.atoms[1].x = r  # 设置距离
        pygcmc.computeSystemEnergy(state)  # 使用默认方法计算能量
        energies_with_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 对比两种方法
    for i, r in enumerate(distances):
        e_hard = energies_hard_cutoff[i]
        e_switch = energies_with_switching[i]
        
        # 小于r_on时，两种方法应该产生相同的能量
        if r < r_on:
            assert e_hard == pytest.approx(e_switch), f"Energy should be the same when r < r_on"
        # 大于r_off时，平滑函数应该使能量为0
        elif r > r_off:
            assert e_switch == pytest.approx(0.0), f"Energy should be 0 when r > r_off"
    
    # 禁用平滑函数，恢复默认状态
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)

def test_with_ewald():
    """
    测试带Ewald求和的平滑函数
    """
    # 创建测试系统
    state = create_test_system()
    
    # 设置原子电荷
    state.atoms[0].charge = 1.0
    state.atoms[1].charge = -1.0
    
    # 设置Ewald参数
    state.info.setTemperature(300.0)  # 300K
    pygcmc.initializeEwaldParameters(state.info.cutoff, state.info.box)
    
    # 设置平滑函数参数
    # 按照 CHARMM 的约定:
    # ctonnb = 1.0 - inner cutoff (r_on)，平滑衰减开始处 (10 Å)
    # ctofnb = 1.2 - outer cutoff (r_off)，势能最终衰减至 0 处 (12 Å)
    r_on = 1.0   # 内截断半径 (10 Å)
    r_off = 1.2  # 外截断半径 (12 Å)
    
    # 设置第二个原子的位置
    distances = [0.9 + i * 0.4/20 for i in range(20)]  # 从0.9到1.3的20个点
    
    # 关闭平滑函数，计算普通Ewald
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)
    
    energies_ewald = []
    for r in distances:
        state.atoms[1].x = r
        pygcmc.computeSystemEnergyEwald(state)  # 使用Ewald方法计算能量
        # 只记录VDW能量部分
        energies_ewald.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 启用平滑函数
    pygcmc.enableSwitchingFunction(state, True, r_on, r_off)
    
    energies_ewald_switching = []
    for r in distances:
        state.atoms[1].x = r
        pygcmc.computeSystemEnergyEwald(state)  # 使用Ewald方法计算能量
        # 只记录VDW能量部分
        energies_ewald_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 对比平滑前后的VDW能量
    for i, r in enumerate(distances):
        e_normal = energies_ewald[i]
        e_switch = energies_ewald_switching[i]
        
        # 小于r_on时，能量应相同
        if r < r_on:
            assert e_normal == pytest.approx(e_switch)
        # 大于r_off时，平滑后VDW能量应为0
        elif r > r_off:
            assert e_switch == pytest.approx(0.0)
    
    # 禁用平滑函数，恢复默认状态
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)

def test_print_switching_values():
    """
    打印平滑函数值和能量数据，用于调试和文档目的
    """
    # 按照 CHARMM 的约定:
    # ctonnb = 1.0 - inner cutoff，平滑衰减开始处 (10 Å)
    # ctofnb = 1.2 - outer cutoff，势能最终衰减至 0 处 (12 Å)
    r_on = 1.0  # 10 Å
    r_off = 1.2  # 12 Å
    
    # 创建测试系统
    state = create_test_system()
    
    # 启用平滑函数
    pygcmc.enableSwitchingFunction(state, True, r_on, r_off)
    
    # 计算平滑函数的值
    key_distances = [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
    switch_values = [pygcmc.calculateSwitchingFunction(state, r) for r in key_distances]
    
    # 打印表头
    print("\nCHARMM Switching Function Values:")
    print(f"{'Distance (nm)':15s} | {'Switch Value':15s}")
    print("-" * 33)
    
    # 打印平滑函数值
    for r, s in zip(key_distances, switch_values):
        print(f"{r:15.3f} | {s:15.4f}")
    
    # 禁用平滑函数，恢复默认状态
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)

def test_compare_energy_with_without_switching():
    """
    详细对比有无switching函数下的能量计算结果
    同时包含direct计算和Ewald计算方法
    """
    # 创建测试系统
    state = create_test_system()
    
    # 设置原子电荷，用于Ewald计算
    state.atoms[0].charge = 1.0
    state.atoms[1].charge = -1.0
    
    # 设置Ewald参数
    state.info.setTemperature(300.0)  # 300K
    pygcmc.initializeEwaldParameters(state.info.cutoff, state.info.box)
    
    # 设置平滑函数参数
    r_on = 1.0   # 内截断半径 (10 Å)
    r_off = 1.2  # 外截断半径 (12 Å)
    
    # 在r_on附近取更多的点，更好展示switching效果
    distances = [0.7 + i * 0.8/50 for i in range(51)]  # 从0.7到1.5的51个点
    
    # 1. Direct计算结果对比
    # 禁用平滑函数
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)
    
    direct_no_switching = []
    for r in distances:
        state.atoms[1].x = r  # 设置距离
        pygcmc.computeSystemEnergy(state)  # 使用direct方法计算能量
        direct_no_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 启用平滑函数
    pygcmc.enableSwitchingFunction(state, True, r_on, r_off)
    
    direct_with_switching = []
    for r in distances:
        state.atoms[1].x = r  # 设置距离
        pygcmc.computeSystemEnergy(state)  # 使用direct方法计算能量
        direct_with_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 2. Ewald计算结果对比
    # 禁用平滑函数
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)
    
    ewald_no_switching_vdw = []
    ewald_no_switching_elec = []
    for r in distances:
        state.atoms[1].x = r  # 设置距离
        pygcmc.computeSystemEnergyEwald(state)  # 使用Ewald方法计算能量
        ewald_no_switching_vdw.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
        ewald_no_switching_elec.append(state.residues[0].energy_elec + state.residues[1].energy_elec)
    
    # 启用平滑函数
    pygcmc.enableSwitchingFunction(state, True, r_on, r_off)
    
    ewald_with_switching_vdw = []
    ewald_with_switching_elec = []
    for r in distances:
        state.atoms[1].x = r  # 设置距离
        pygcmc.computeSystemEnergyEwald(state)  # 使用Ewald方法计算能量
        ewald_with_switching_vdw.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
        ewald_with_switching_elec.append(state.residues[0].energy_elec + state.residues[1].energy_elec)
    
    # 打印结果并验证
    print("\n能量计算对比结果 (带switching vs 不带switching):")
    print(f"{'距离(nm)':10s} | {'Direct无切换':14s} | {'Direct有切换':14s} | {'Ewald VDW无切换':16s} | {'Ewald VDW有切换':16s}")
    print("-" * 80)
    
    # 打印部分关键结果点
    key_indices = [0, 10, 20, 25, 30, 35, 40, 45, 50]  # 选择几个关键点展示
    
    for idx in key_indices:
        if idx < len(distances):
            r = distances[idx]
            # 在r_off处计算平滑函数值便于验证
            switch_value = pygcmc.calculateSwitchingFunction(state, r) if r_on <= r <= r_off else (1.0 if r < r_on else 0.0)
            
            print(f"{r:10.3f} | {direct_no_switching[idx]:14.6f} | {direct_with_switching[idx]:14.6f} | "
                  f"{ewald_no_switching_vdw[idx]:16.6f} | {ewald_with_switching_vdw[idx]:16.6f}")
    
    # 进行验证
    for i, r in enumerate(distances):
        # 1. 距离小于r_on时，能量应该相同
        if r < r_on:
            assert direct_no_switching[i] == pytest.approx(direct_with_switching[i]), \
                  f"Direct energy should be the same when r < r_on, at r={r}"
            assert ewald_no_switching_vdw[i] == pytest.approx(ewald_with_switching_vdw[i]), \
                  f"Ewald VDW energy should be the same when r < r_on, at r={r}"
                  
        # 2. 距离大于r_off时，使用switching的能量应为0
        elif r > r_off:
            assert direct_with_switching[i] == pytest.approx(0.0), \
                  f"Direct energy with switching should be 0 when r > r_off, at r={r}"
            assert ewald_with_switching_vdw[i] == pytest.approx(0.0), \
                  f"Ewald VDW energy with switching should be 0 when r > r_off, at r={r}"
                  
        # 3. 在r_on和r_off之间，检查switching是否正确应用
        elif r_on <= r <= r_off:
            # 计算预期的switching值
            switch_value = pygcmc.calculateSwitchingFunction(state, r)
            
            # 检查direct能量是否按switching缩放
            assert direct_with_switching[i] == pytest.approx(direct_no_switching[i] * switch_value, abs=1e-5), \
                  f"Direct energy with switching not correctly scaled at r={r}"
                  
            # 检查Ewald VDW能量是否按switching缩放
            assert ewald_with_switching_vdw[i] == pytest.approx(ewald_no_switching_vdw[i] * switch_value, abs=1e-5), \
                  f"Ewald VDW energy with switching not correctly scaled at r={r}"
    
    # 电荷能量测试 - 只检查关键点
    if ewald_no_switching_elec[0] != 0:  # 确保有电荷能量
        print("\n静电能量对比结果 (Ewald):")
        print(f"{'距离(nm)':10s} | {'Ewald静电无切换':18s} | {'Ewald静电有切换':18s}")
        print("-" * 60)
        
        for idx in key_indices:
            if idx < len(distances):
                r = distances[idx]
                print(f"{r:10.3f} | {ewald_no_switching_elec[idx]:18.6f} | {ewald_with_switching_elec[idx]:18.6f}")
    
    # 禁用平滑函数，恢复默认状态
    pygcmc.enableSwitchingFunction(state, False, r_on, r_off)

if __name__ == "__main__":
    # 运行测试函数
    test_print_switching_values()
    # 运行新增的对比测试
    test_compare_energy_with_without_switching() 