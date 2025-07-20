#!/usr/bin/env python3
"""
测试FBP的物理正确性
"""

import numpy as np

def test_physics():
    """测试基本物理"""
    print("基本物理测试")
    print("="*60)
    
    # 情况1：正电荷在右，负电荷在左
    print("\n情况1：正电荷在(1,0,0)，计算在(0,0,0)处的电场和负电荷受力")
    
    # 电场计算
    q_source = 1.0  # 正电荷
    r_source = np.array([1.0, 0.0, 0.0])
    r_field = np.array([0.0, 0.0, 0.0])
    
    # 从源到场点的矢量
    r_vec = r_field - r_source  # (-1, 0, 0)
    r_mag = np.linalg.norm(r_vec)  # 1.0
    r_hat = r_vec / r_mag  # (-1, 0, 0)
    
    # 电场：E = k*q/r^2 * r_hat（从正电荷指向外）
    k = 138.935  # kJ/(mol·nm·e^2)
    E = k * q_source / r_mag**2 * r_hat
    print(f"  r_vec = {r_vec}")
    print(f"  r_hat = {r_hat}")
    print(f"  电场 E = {E} kJ/(mol·nm·e)")
    
    # 负电荷受力
    q_test = -1.0
    F = q_test * E
    print(f"  负电荷受力 F = q*E = {F} kJ/(mol·nm)")
    print(f"  力的方向：{'向右（被吸引）' if F[0] > 0 else '向左（被排斥）'}")
    
    # 弹簧平衡
    k_spring = 138935.0  # kJ/(mol·nm^2)
    displacement = F / k_spring
    print(f"  平衡位移 = F/k = {displacement} nm = {displacement*1000} pm")
    
    print("\n物理检查：")
    print(f"  ✓ 负电荷应该被正电荷吸引，向右移动")
    print(f"  {'✓' if F[0] > 0 else '✗'} 实际力方向{'正确' if F[0] > 0 else '错误'}")
    print(f"  {'✓' if displacement[0] > 0 else '✗'} 实际位移方向{'正确' if displacement[0] > 0 else '错误'}")
    
    # 情况2：FBP算法的计算方式
    print("\n\n情况2：模拟FBP算法的计算")
    
    # 按照calculateFixedElectricField的方式
    # double dx = drude.x - atom.x;  // r_field - r_source
    dx = 0.0 - 1.0  # -1.0
    r = abs(dx)  # 1.0
    factor = q_source / (r**3)  # 1.0
    E_fbp_x = factor * dx  # -1.0
    
    print(f"  dx = drude.x - atom.x = {dx}")
    print(f"  factor = charge / r^3 = {factor}")
    print(f"  E_x = factor * dx = {E_fbp_x}")
    
    # 转换为力
    F_fbp = q_test * k * E_fbp_x
    print(f"  F = q * k * E = {q_test} * {k} * {E_fbp_x} = {F_fbp}")
    
    # 计算位移
    displacement_fbp = F_fbp / k_spring
    print(f"  位移 = {displacement_fbp} nm = {displacement_fbp*1000} pm")
    
    print("\n\n结论：")
    print("FBP中的电场计算可能缺少了单位转换因子k=138.935")

if __name__ == "__main__":
    test_physics()