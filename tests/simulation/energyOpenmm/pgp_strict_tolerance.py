"""
测试 PGP 在严格容错标准下的表现

重现原始 PGP 测试场景，但使用更严格的容错标准
"""

import sys
import os
import random
import statistics

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import (setPGPParameters, initializePMEParameters, 
                    precomputeGridPotential, computeSystemEnergyPGP,
                    calculateMoleculeEnergy, computeMovementEnergyPME,
                    computeSystemEnergyPME, computeSystemEnergyEwald)


def test_pgp_with_strict_tolerance():
    """测试 PGP 在不同容错标准下的表现"""
    
    print("\n" + "="*80)
    print("PGP 严格容错测试")
    print("="*80)
    
    # 创建一个简单的测试系统
    box_size = 5.0  # nm
    cutoff = 1.8    # nm
    
    # Fixed molecules
    fixed_positions = [
        [1.0, 1.0, 2.5],
        [4.0, 1.0, 2.5],
        [4.0, 4.0, 2.5],
        [1.0, 4.0, 2.5]
    ]
    fixed_charges = [1.0, -1.0, 1.0, -1.0]
    
    # Moving molecule (initial position)
    moving_position_initial = [2.5, 2.5, 2.5]
    moving_charge = 0.5
    
    # 创建系统
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]    # 无 LJ，只测试静电
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Add atoms
    all_positions_initial = fixed_positions + [moving_position_initial]
    all_charges = fixed_charges + [moving_charge]
    
    atoms = []
    for i in range(5):
        atom = MCAtom()
        atom.x, atom.y, atom.z = all_positions_initial[i]
        atom.charge = all_charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    
    # 添加残基
    residues = []
    for i in range(5):
        res = MCResidue()
        res.active = True
        res.fixed = (i < 4)  # 前 4 个固定
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 5
    
    # 初始化参数
    alpha = 2.5
    mesh_size = [32, 32, 32]
    
    initializePMEParameters(cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, cutoff, mesh_size, 4, 1e-6)
    
    # 初始化 Ewald 参数
    pygcmc.initializeEwaldParameters(cutoff, state.info.box, alpha)
    
    # 预计算网格
    precomputeGridPotential(state, fixed_only=True)
    
    # 设置移动残基
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 4
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    print("\n测试多个随机移动...")
    print("-" * 60)
    
    errors_pgp_pme = []
    errors_pgp_ewald = []
    n_moves = 10
    
    for move_idx in range(n_moves):
        # 计算初始能量
        initial_pgp = calculateMoleculeEnergy(state)
        initial_pme_result = computeMovementEnergyPME(state)
        initial_pme_recip = initial_pme_result[2].get('reciprocal', 0.0)
        initial_ewald_result = computeSystemEnergyEwald(state)
        initial_ewald_recip = initial_ewald_result[2].get('reciprocal', 0.0)
        
        # Random movement (small displacement)
        dx = random.uniform(-0.3, 0.3)
        dy = random.uniform(-0.3, 0.3)
        dz = random.uniform(-0.3, 0.3)
        
        # 移动原子
        moving_atom = state.atoms[4]
        moving_atom.x = (moving_position_initial[0] + dx) % box_size
        moving_atom.y = (moving_position_initial[1] + dy) % box_size
        moving_atom.z = (moving_position_initial[2] + dz) % box_size
        
        # 计算移动后能量
        moved_pgp = calculateMoleculeEnergy(state)
        moved_pme_result = computeMovementEnergyPME(state)
        moved_pme_recip = moved_pme_result[2].get('reciprocal', 0.0)
        moved_ewald_result = computeSystemEnergyEwald(state)
        moved_ewald_recip = moved_ewald_result[2].get('reciprocal', 0.0)
        
        # 计算能量变化
        delta_pgp = moved_pgp - initial_pgp
        delta_pme = moved_pme_recip - initial_pme_recip
        delta_ewald = moved_ewald_recip - initial_ewald_recip
        
        print(f"\nMove {move_idx + 1}: ({dx:.3f}, {dy:.3f}, {dz:.3f})")
        print(f"ΔE_PGP:   {delta_pgp:10.6f} kJ/mol")
        print(f"ΔE_PME:   {delta_pme:10.6f} kJ/mol")
        print(f"ΔE_Ewald: {delta_ewald:10.6f} kJ/mol")
        
        # 计算相对误差
        if abs(delta_pme) > 1e-6:
            error_pgp_pme = abs((delta_pgp - delta_pme) / delta_pme)
            errors_pgp_pme.append(error_pgp_pme)
            print(f"PGP vs PME 误差: {error_pgp_pme*100:.2f}%")
        
        if abs(delta_ewald) > 1e-6:
            error_pgp_ewald = abs((delta_pgp - delta_ewald) / delta_ewald)
            errors_pgp_ewald.append(error_pgp_ewald)
            print(f"PGP vs Ewald 误差: {error_pgp_ewald*100:.2f}%")
        
        # 恢复原始位置
        moving_atom.x = moving_position_initial[0]
        moving_atom.y = moving_position_initial[1]
        moving_atom.z = moving_position_initial[2]
    
    # 分析结果
    print("\n" + "="*80)
    print("误差分析")
    print("="*80)
    
    if errors_pgp_pme:
        avg_error = statistics.mean(errors_pgp_pme) * 100
        max_error = max(errors_pgp_pme) * 100
        min_error = min(errors_pgp_pme) * 100
        
        print(f"\nPGP vs PME 误差统计:")
        print(f"平均误差: {avg_error:.2f}%")
        print(f"最大误差: {max_error:.2f}%")
        print(f"最小误差: {min_error:.2f}%")
        
        # 测试不同的容错标准
        tolerances = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
        print(f"\n不同容错标准下的测试结果:")
        print("-" * 40)
        
        for tol in tolerances:
            passed = avg_error/100 < tol
            status = "✓ Pass" if passed else "✗ Fail"
            print(f"容错 {tol*100:5.1f}%: {status}")
        
        # 检查具体有多少测试会失败
        print(f"\n单个测试的通过率:")
        for tol in tolerances:
            n_passed = sum(1 for e in errors_pgp_pme if e < tol)
            pass_rate = n_passed / len(errors_pgp_pme) * 100
            print(f"容错 {tol*100:5.1f}%: {n_passed}/{len(errors_pgp_pme)} ({pass_rate:.1f}%)")
    
    # 测试系统总能量
    print("\n" + "="*80)
    print("系统总能量测试")
    print("="*80)
    
    # 计算系统总能量
    computeSystemEnergyPGP(state)
    pgp_total = state.ewald_energy.get('total', 0.0)
    pgp_recip = state.ewald_energy.get('reciprocal', 0.0)
    
    computeSystemEnergyPME(state)
    pme_total = state.ewald_energy.get('total', 0.0)
    pme_recip = state.ewald_energy.get('reciprocal', 0.0)
    
    print(f"\nPME 倒空间: {pme_recip:.2f} kJ/mol")
    print(f"PGP 倒空间: {pgp_recip:.2f} kJ/mol")
    print(f"比例 PGP/PME: {pgp_recip/pme_recip if pme_recip != 0 else 0:.3f}")
    
    if abs(pgp_recip/pme_recip - 2.0) < 0.01:
        print("\n⚠️  发现 PGP 倒空间能量是 PME 的 2 倍！")
    
    print("\n" + "="*80)
    print("结论")
    print("="*80)
    
    if errors_pgp_pme and avg_error > 10:
        print("✗ PGP 在严格容错标准下无法通过测试")
        print(f"  平均误差 {avg_error:.1f}% 远超合理范围")
    else:
        print("? 需要更多测试数据")


if __name__ == "__main__":
    test_pgp_with_strict_tolerance()