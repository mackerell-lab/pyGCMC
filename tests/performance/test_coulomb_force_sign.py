#!/usr/bin/env python3
"""
测试库仑力的符号约定
"""

import numpy as np

def analyze_coulomb_force_calculation():
    """分析库仑力计算"""
    print("库仑力计算分析")
    print("="*60)
    
    print("\n查看DrudeForce.cpp第526-531行：")
    print("```cpp")
    print("// Force")
    print("double forceMag = ONE_4PI_EPS0 * chargeProduct * invR2;")
    print("Vec3 f = delta * (forceMag * invR);")
    print("")
    print("forces[i] += f;")
    print("forces[j] -= f;")
    print("```")
    
    print("\n符号分析：")
    print("-"*40)
    
    # 例子：i是Drude(负电荷)，j是External(正电荷)
    print("例子：i=Drude(-1), j=External(+1)")
    print("位置：Drude在(0,0,0), External在(1,0,0)")
    
    # delta = pos2 - pos1 = pos_j - pos_i
    delta_x = 1.0 - 0.0  # = 1.0
    print(f"\ndelta = pos_j - pos_i = (1,0,0)")
    print(f"r = |delta| = 1.0")
    
    # 力的大小
    k = 138.935
    q_i = -1.0
    q_j = 1.0
    chargeProduct = q_i * q_j  # = -1.0
    forceMag = k * chargeProduct / 1.0**2  # = -138.935
    
    print(f"\nchargeProduct = {chargeProduct}")
    print(f"forceMag = k * chargeProduct / r² = {forceMag}")
    
    # 力的方向
    # f = delta * (forceMag / r) = (1,0,0) * (-138.935 / 1.0)
    f_x = delta_x * (forceMag / 1.0)  # = 1.0 * (-138.935) = -138.935
    
    print(f"\nf = delta * (forceMag / r)")
    print(f"f_x = {delta_x} * ({forceMag} / 1.0) = {f_x}")
    
    print(f"\nforces[Drude] += f  →  forces[Drude] += {f_x}")
    print(f"forces[External] -= f  →  forces[External] -= {f_x}")
    
    print("\n物理解释：")
    print(f"- Drude受力: {f_x} kJ/(mol·nm) (向左)")
    print(f"- 这是错误的！负电荷应该被正电荷吸引，力应该向右(正)")
    
    print("\n\n问题诊断：")
    print("="*60)
    print("库仑力的符号约定有问题！")
    print("")
    print("正确的计算应该是：")
    print("- 从i指向j的矢量：r_ij = r_j - r_i")
    print("- i受到的力：F_i = k*q_i*q_j/r² * r_hat_ij")
    print("- 对于吸引力(q_i*q_j < 0)，力应该沿着r_ij方向")
    print("")
    print("但代码中：")
    print("- f = delta * (k*q_i*q_j/r³)")
    print("- 当q_i*q_j < 0时，f与delta反向")
    print("- 这导致吸引力变成了排斥力！")

def show_correct_implementation():
    """展示正确的实现"""
    print("\n\n正确的实现应该是：")
    print("="*60)
    
    print("```cpp")
    print("// 正确的库仑力计算")
    print("Vec3 r_ij = pos_j - pos_i;  // 从i指向j")
    print("double r = r_ij.norm();")
    print("Vec3 r_hat = r_ij / r;      // 单位矢量")
    print("")
    print("// 库仑定律：F_i = k*q_i*q_j/r² * r_hat")
    print("double F_mag = k * q_i * q_j / (r * r);")
    print("Vec3 F_i = r_hat * F_mag;")
    print("")
    print("// 如果q_i*q_j < 0 (异号电荷):")
    print("// - F_mag < 0")
    print("// - F_i 与 r_hat 反向")
    print("// - 即力指向对方（吸引）")
    print("```")
    
    print("\n但SCF代码使用了简化形式，可能导致符号错误")

def compare_scf_fbp_approach():
    """比较SCF和FBP的方法"""
    print("\n\nSCF vs FBP方法对比")
    print("="*60)
    
    print("SCF方法：")
    print("1. 计算所有力（弹簧+库仑）")
    print("2. 更新位置：r_new = r_old + damping * F_total / k")
    print("3. 问题：库仑力符号错误，导致错误的更新方向")
    print("")
    
    print("FBP方法：")
    print("1. 计算电场（而非力）")
    print("2. 求解平衡方程：F_spring + q*E = 0")
    print("3. 直接得到位移：d = q*E / k")
    print("4. 避免了力的符号混淆")
    
    print("\n\n结论：")
    print("-"*40)
    print("SCF的bug在于calculateCoulombEnergy中的符号约定！")
    print("当计算异号电荷间的吸引力时，力的方向错误。")
    print("这解释了为什么SCF会向错误的方向收敛。")

if __name__ == "__main__":
    analyze_coulomb_force_calculation()
    show_correct_implementation()
    compare_scf_fbp_approach()