#!/usr/bin/env python3
"""
调试FBP中的电场和力计算
"""

import numpy as np

def debug_fbp_calculation():
    """调试FBP的计算逻辑"""
    print("FBP算法调试")
    print("="*60)
    
    # 系统设置
    # Drude在(0,0,0)，外部正电荷在(1,0,0)
    drude_pos = np.array([0.0, 0.0, 0.0])
    external_pos = np.array([1.0, 0.0, 0.0])
    
    q_drude = -1.0
    q_external = 1.0
    k_elec = 138.935
    
    print("系统配置:")
    print(f"  Drude位置: {drude_pos}")
    print(f"  Drude电荷: {q_drude}")
    print(f"  外部电荷位置: {external_pos}")
    print(f"  外部电荷: {q_external}")
    
    # 按照DrudeForceBalance.cpp中calculateFixedElectricField的方式计算
    print("\n按照FBP代码计算电场:")
    print("-"*40)
    
    # double dx = drude.x - atom.x;  // r_field - r_source
    dx = drude_pos[0] - external_pos[0]  # 0 - 1 = -1
    dy = drude_pos[1] - external_pos[1]  # 0
    dz = drude_pos[2] - external_pos[2]  # 0
    
    print(f"  dx = drude.x - atom.x = {dx}")
    
    r2 = dx*dx + dy*dy + dz*dz  # 1.0
    r = np.sqrt(r2)  # 1.0
    
    # 原代码：double factor = 138.935 * atom.charge / (r2 * r);
    factor = k_elec * q_external / (r2 * r)  # 138.935 * 1.0 / 1.0 = 138.935
    
    print(f"  r = {r}")
    print(f"  factor = k*q/(r^3) = {factor}")
    
    # 电场分量
    E_x = factor * dx  # 138.935 * (-1) = -138.935
    E_y = factor * dy  # 0
    E_z = factor * dz  # 0
    
    print(f"  电场 E = ({E_x}, {E_y}, {E_z})")
    print(f"  |E| = {np.sqrt(E_x**2 + E_y**2 + E_z**2)}")
    
    # 力计算（按照calculateDrudeForces）
    F_x = q_drude * E_x  # -1 * (-138.935) = 138.935
    F_y = q_drude * E_y  # 0
    F_z = q_drude * E_z  # 0
    
    print(f"\n力计算:")
    print(f"  F = q*E = {q_drude} * ({E_x}, {E_y}, {E_z})")
    print(f"  F = ({F_x}, {F_y}, {F_z}) kJ/(mol·nm)")
    
    # 物理分析
    print("\n物理分析:")
    print("-"*40)
    print(f"  外部正电荷在右侧(+x方向)")
    print(f"  Drude负电荷应该被吸引，向右移动")
    print(f"  因此力应该指向+x方向（正值）")
    print(f"  计算得到的力: F_x = {F_x} ✓ 正确")
    
    # 位移计算
    k_spring = 138935.0
    displacement = F_x / k_spring
    print(f"\n位移计算:")
    print(f"  位移 = F/k = {F_x}/{k_spring} = {displacement:.6f} nm")
    print(f"  位移 = {displacement*1000:.3f} pm")
    
    # 等等！让我们仔细想想
    print("\n\n深入分析:")
    print("="*60)
    print("FBP计算出F_x = +138.935，这意味着力指向+x（向右）")
    print("这对于负电荷被正电荷吸引是正确的")
    print("位移 = F/k = +0.001 nm也是正确的（向右）")
    print("")
    print("但测试显示FBP给出+1.0 pm的位移，而SCF给出-0.998 pm")
    print("这说明问题可能在于坐标系或者力的定义")
    
    # 重新分析
    print("\n\n重新分析SCF vs FBP:")
    print("-"*60)
    print("SCF: 位移 = -0.998 pm")
    print("  这意味着Drude向左移动（x坐标减小）")
    print("  但这物理上是错误的！负电荷应该被右侧的正电荷吸引")
    print("")
    print("FBP: 位移 = +1.000 pm")  
    print("  这意味着Drude向右移动（x坐标增加）")
    print("  这物理上是正确的！负电荷被右侧的正电荷吸引")
    print("")
    print("结论：可能是SCF的结果解释有误，或者测试代码有问题")

if __name__ == "__main__":
    debug_fbp_calculation()