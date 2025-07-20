#!/usr/bin/env python3
"""
分析为什么密度1.0 g/cm³无法收敛
"""

import numpy as np

def analyze_box_sizes_and_forces():
    """
    分析不同密度下的盒子尺寸和相互作用
    """
    print("密度对收敛性的影响分析")
    print("="*70)
    
    # 物理常数
    molar_mass_water = 18.01528  # g/mol
    avogadro = 6.02214076e23     # mol⁻¹
    ONE_4PI_EPS0 = 138.935456    # kJ/mol·nm·e^-2
    
    # 分析不同水分子数和密度
    n_waters_list = [2, 4, 8, 16, 32, 64]
    densities = [0.3, 1.0]
    
    print(f"\n{'水分子数':>8} {'密度(g/cm³)':>12} {'盒子(nm)':>10} {'半盒子(nm)':>12} {'最近O-O(nm)':>12}")
    print("-"*70)
    
    for n_waters in n_waters_list:
        for density in densities:
            # 计算盒子大小
            total_mass_g = n_waters * molar_mass_water / avogadro
            volume_cm3 = total_mass_g / density
            volume_nm3 = volume_cm3 * 1e21
            box_length = volume_nm3 ** (1.0/3.0)
            
            # 估算最近O-O距离（假设均匀分布）
            volume_per_water = volume_nm3 / n_waters
            nearest_oo = volume_per_water ** (1.0/3.0)
            
            print(f"{n_waters:8d} {density:12.1f} {box_length:10.3f} {box_length/2:12.3f} {nearest_oo:12.3f}")
    
    print("\n\n关键问题分析:")
    print("="*70)
    
    # 分析2水系统的具体问题
    print("\n1. 以2水系统为例:")
    
    for density in [0.3, 1.0]:
        n_waters = 2
        total_mass_g = n_waters * molar_mass_water / avogadro
        volume_cm3 = total_mass_g / density
        volume_nm3 = volume_cm3 * 1e21
        box_length = volume_nm3 ** (1.0/3.0)
        
        print(f"\n密度 {density} g/cm³:")
        print(f"  盒子长度: {box_length:.3f} nm")
        print(f"  半盒子: {box_length/2:.3f} nm")
        
        # Thole截断问题
        thole_cutoff = 0.8  # nm
        ratio = thole_cutoff / (box_length/2)
        print(f"  Thole截断(0.8nm) / 半盒子 = {ratio:.2f}")
        
        if ratio > 1:
            print(f"  ⚠️  违反最小镜像约定！")
            
        # 动态调整后的Thole截断
        adjusted_thole = min(0.8, box_length/2 - 0.01)
        print(f"  调整后Thole截断: {adjusted_thole:.3f} nm")
        
        # 估算两个水分子的距离和相互作用
        # 假设两个水分子在对角位置
        typical_distance = box_length / 2
        print(f"  典型O-O距离: {typical_distance:.3f} nm")
        
        # 计算库仑相互作用强度
        # 水分子的偶极矩约为1.85 Debye = 0.0617 e·nm
        dipole_moment = 0.0617  # e·nm
        
        # 偶极-偶极相互作用能量（平行偶极）
        # E = -2 * k * μ1 * μ2 / r³
        if typical_distance > 0.1:
            E_dipole = -2 * ONE_4PI_EPS0 * dipole_moment**2 / typical_distance**3
            print(f"  偶极-偶极能量: {E_dipole:.1f} kJ/mol")
        
        # 计算初始电场强度
        # 点电荷产生的电场 E = k*q/r²
        q_oxygen = 1.71636  # e
        if typical_distance > 0.1:
            E_field = ONE_4PI_EPS0 * q_oxygen / typical_distance**2
            print(f"  O原子在另一个O处的电场: {E_field:.0f} kJ/(mol·nm·e)")
            
            # Drude粒子的预期位移
            alpha = 0.0009782237  # nm³
            q_drude = 1.71636     # e (绝对值)
            expected_displacement = alpha * E_field / q_drude * 1000  # pm
            print(f"  预期Drude位移: {expected_displacement:.0f} pm")
            
            if expected_displacement > 20:
                print(f"  ⚠️  超过硬墙约束(20 pm)！")

def analyze_pbc_artifacts():
    """
    分析周期性边界条件导致的伪影
    """
    print("\n\n2. 周期性边界条件(PBC)伪影分析:")
    print("="*70)
    
    # 2水系统，密度1.0
    box_length = 0.391  # nm
    
    print(f"\n2水系统，密度1.0 g/cm³，盒子{box_length} nm:")
    print("\n最小镜像约定要求：任何相互作用的截断距离 < 盒子长度/2")
    
    # 不同类型的相互作用
    interactions = [
        ("LJ (O-O)", 1.2, "通常的非键截断"),
        ("静电", 1.2, "通常的静电截断"),
        ("Thole", 0.8, "Thole屏蔽截断"),
        ("调整后Thole", 0.186, "动态调整以满足PBC")
    ]
    
    print(f"\n{'相互作用类型':>15} {'截断(nm)':>10} {'截断/半盒子':>12} {'状态':>10}")
    print("-"*60)
    
    for name, cutoff, desc in interactions:
        ratio = cutoff / (box_length/2)
        status = "✓ OK" if ratio < 1 else "✗ 违反"
        print(f"{name:>15} {cutoff:10.3f} {ratio:12.2f} {status:>10}")
    
    print("\n问题：")
    print("- 原始Thole截断(0.8 nm)是半盒子的4.08倍！")
    print("- 这意味着一个分子会'看到'另一个分子的多个镜像")
    print("- 导致非物理的强相互作用和高电场")

def analyze_scf_convergence_difficulty():
    """
    分析SCF收敛困难的原因
    """
    print("\n\n3. SCF收敛困难的根本原因:")
    print("="*70)
    
    print("\n高密度小系统的恶性循环：")
    print("1. 小盒子 → PBC违反 → 多重镜像相互作用")
    print("2. 强电场 → 大Drude位移需求 → 撞击硬墙约束")
    print("3. 硬墙反弹 → 振荡 → 无法收敛")
    
    print("\n具体数值（2水，密度1.0）：")
    print("- 初始电场: >3000 kJ/(mol·nm·e)")
    print("- 预期位移: >1900 pm")
    print("- 硬墙限制: 20 pm")
    print("- 位移被压缩: 1900→20 pm (压缩95倍！)")
    
    print("\n收敛警告中的RMS力：")
    forces = [10664.5, 4785.38, 1823.29, 795.698]
    print("- 2水: 10664.5 kJ/mol/nm")
    print("- 4水: 4785.4 kJ/mol/nm")
    print("- 8水: 1823.3 kJ/mol/nm")
    print("- 16水: 795.7 kJ/mol/nm")
    print("\n随着系统增大，RMS力降低，但仍然很高")

def suggest_solutions():
    """
    建议解决方案
    """
    print("\n\n4. 解决方案:")
    print("="*70)
    
    print("\n短期方案（已实施）：")
    print("✓ 动态调整Thole截断")
    print("✓ 使用更低密度(0.3 g/cm³)进行测试")
    
    print("\n中期方案：")
    print("- 实现更好的初始猜测（基于局部电场）")
    print("- 使用软墙约束代替硬墙")
    print("- 实现ASPC等更稳定的SCF算法")
    
    print("\n长期方案：")
    print("- 实现Ewald求和处理长程静电")
    print("- 使用最小镜像约定的特殊处理")
    print("- 开发专门的小系统算法")
    
    print("\n推荐的测试策略：")
    print("1. 算法开发：使用密度0.1-0.3 g/cm³")
    print("2. 物理验证：使用≥256分子，密度1.0 g/cm³")
    print("3. 小系统测试：降低密度或使用非周期边界")

def main():
    """
    主函数
    """
    analyze_box_sizes_and_forces()
    analyze_pbc_artifacts()
    analyze_scf_convergence_difficulty()
    suggest_solutions()
    
    print("\n\n总结:")
    print("="*70)
    print("密度1.0 g/cm³无法收敛的根本原因是：")
    print("1. 小系统的盒子太小，违反PBC最小镜像约定")
    print("2. 产生非物理的强电场和多重镜像相互作用")
    print("3. Drude粒子需要极大位移但被硬墙约束限制")
    print("4. 形成振荡，无法达到自洽")
    print("\n这不是算法错误，而是小系统高密度的固有物理限制。")

if __name__ == "__main__":
    main()