#!/usr/bin/env python3
"""
测试电场方向约定
"""

import numpy as np

def test_electric_field_convention():
    """测试电场方向的约定"""
    print("电场方向约定测试")
    print("="*60)
    
    # 设置
    q_source = 1.0  # 正电荷
    r_source = np.array([1.0, 0.0, 0.0])  # 在x=1处
    r_field = np.array([0.0, 0.0, 0.0])   # 在原点
    
    print("配置:")
    print(f"  源电荷: q = +{q_source} at {r_source}")
    print(f"  场点: {r_field}")
    
    # 方法1：标准物理学定义
    print("\n方法1：标准物理学定义")
    print("-"*40)
    # 从源指向场点的矢量
    r_vec = r_field - r_source  # (-1, 0, 0)
    r_mag = np.linalg.norm(r_vec)
    r_hat = r_vec / r_mag
    
    # 电场
    k = 138.935
    E_standard = k * q_source / r_mag**2 * r_hat
    
    print(f"  从源到场点: r = {r_vec}")
    print(f"  单位矢量: r_hat = {r_hat}")
    print(f"  电场: E = {E_standard} kJ/(mol·nm·e)")
    print(f"  解释：正电荷在原点产生向左的电场")
    
    # 测试电荷受力
    q_test = -1.0
    F_standard = q_test * E_standard
    print(f"\n  测试电荷 q = {q_test}:")
    print(f"  受力 F = q*E = {F_standard}")
    print(f"  解释：负电荷在向左的电场中受向右的力（被吸引）")
    
    # 方法2：FBP代码中的实现
    print("\n\n方法2：FBP代码实现")
    print("-"*40)
    # dx = drude.x - atom.x
    dx = r_field[0] - r_source[0]  # 0 - 1 = -1
    r = abs(dx)
    factor = k * q_source / (r**3)
    E_fbp_x = factor * dx
    
    print(f"  dx = field.x - source.x = {dx}")
    print(f"  factor = k*q/r³ = {factor}")
    print(f"  E_x = factor * dx = {E_fbp_x}")
    
    F_fbp = q_test * E_fbp_x
    print(f"\n  测试电荷 q = {q_test}:")
    print(f"  受力 F = q*E = {F_fbp}")
    
    # 比较
    print("\n\n比较:")
    print("="*40)
    print(f"标准方法: E = {E_standard[0]:.1f}, F = {F_standard[0]:.1f}")
    print(f"FBP方法:  E = {E_fbp_x:.1f}, F = {F_fbp:.1f}")
    print(f"结果: {'一致' if abs(E_standard[0] - E_fbp_x) < 0.1 else '不一致'}")
    
    # 关键洞察
    print("\n关键洞察:")
    print("-"*40)
    print("两种方法给出相同的结果！")
    print("FBP的电场计算是正确的。")
    print("")
    print("那么为什么FBP和SCF给出相反的位移？")
    print("")
    print("可能的原因：")
    print("1. 初始猜测的差异")
    print("2. 迭代更新的方式不同")
    print("3. 力平衡的定义不同")
    print("4. Parent电荷的影响（虽然应该被排除）")
    
    # 深入分析
    print("\n\n深入分析：为什么结果相反？")
    print("="*60)
    
    print("\n情况1：Parent不带电")
    print("- 理论：Drude应向左移动（远离右侧正电荷）")
    print("- SCF：向左 ✓ 符合物理")
    print("- FBP：向右 ✗ 不符合物理")
    
    print("\n情况2：Parent带正电")
    print("- Parent和Drude相互吸引（但应被排除）")
    print("- 如果没有正确排除，Drude会被Parent强烈吸引")
    print("- 这可能解释为什么SCF仍然向左（可能包含了Parent的影响）")
    
    print("\n猜测：")
    print("1. SCF可能在计算中包含了Parent-Drude相互作用")
    print("2. FBP的力平衡逻辑可能有符号错误")
    print("3. 两种算法优化的目标函数可能不同")

if __name__ == "__main__":
    test_electric_field_convention()