#!/usr/bin/env python3
"""
重新分析SCF的力方向问题
"""

import numpy as np

def analyze_physics_conventions():
    """分析物理约定"""
    print("物理约定分析")
    print("="*60)
    
    print("\n标准库仑定律：")
    print("F_ij = k * q_i * q_j / r² * r_hat_ij")
    print("其中 r_hat_ij 是从i指向j的单位矢量")
    print("")
    print("如果我们定义 delta = r_j - r_i（从i指向j），那么：")
    print("r_hat = delta / |delta|")
    print("")
    print("所以：F_ij = k * q_i * q_j / |delta|³ * delta")
    
    print("\n\n例子验证：")
    print("-"*40)
    
    # 情况1：两个正电荷
    print("\n情况1：q_i = +1, q_j = +1")
    print("位置：i在原点，j在(1,0,0)")
    q_i, q_j = 1.0, 1.0
    delta = np.array([1.0, 0.0, 0.0])
    force_factor = q_i * q_j / np.linalg.norm(delta)**3
    F = force_factor * delta
    print(f"F_i = {force_factor} * {delta} = {F}")
    print("解释：i受到向右的力（+x方向），被j排斥 ✓")
    
    # 情况2：异号电荷
    print("\n情况2：q_i = -1, q_j = +1")
    q_i, q_j = -1.0, 1.0
    force_factor = q_i * q_j / np.linalg.norm(delta)**3
    F = force_factor * delta
    print(f"F_i = {force_factor} * {delta} = {F}")
    print("解释：i受到向左的力（-x方向），被j吸引？✗")
    print("错误！负电荷应该被正电荷吸引，力应该向右！")
    
    print("\n\n问题分析：")
    print("="*60)
    print("等等！我意识到了问题所在！")
    print("")
    print("在物理学中，库仑力的完整表达式是：")
    print("F_12 = k * q_1 * q_2 / r² * r_hat_12")
    print("")
    print("但这里的r_hat_12是有歧义的！通常有两种约定：")
    print("")
    print("约定1（物理教科书）：")
    print("- r_12 = r_2 - r_1（从1指向2）")
    print("- F_12 = 作用在粒子1上的力")
    print("- 对于q1<0, q2>0：F_12沿着r_12方向（吸引）")
    print("")
    print("约定2（某些代码实现）：")
    print("- 可能使用不同的符号约定")

def check_our_implementation():
    """检查我们的实现"""
    print("\n\n我们的实现分析")
    print("="*60)
    
    print("在DrudeForce.cpp中：")
    print("```cpp")
    print("Vec3 delta = pos2 - pos1;  // 从i指向j")
    print("double forceMag = ONE_4PI_EPS0 * chargeProduct * invR2;")
    print("Vec3 f = delta * (forceMag * invR);")
    print("forces[i] += f;")
    print("```")
    print("")
    print("让我们用数值验证：")
    
    # Drude在原点，External在(1,0,0)
    pos_drude = np.array([0.0, 0.0, 0.0])
    pos_external = np.array([1.0, 0.0, 0.0])
    q_drude = -1.0
    q_external = 1.0
    k = 138.935
    
    delta = pos_external - pos_drude  # (1,0,0)
    r = np.linalg.norm(delta)
    chargeProduct = q_drude * q_external  # -1
    forceMag = k * chargeProduct / r**2  # -138.935
    f = delta * (forceMag / r)  # (1,0,0) * (-138.935) = (-138.935,0,0)
    
    print(f"\ndelta = {delta}")
    print(f"chargeProduct = {chargeProduct}")
    print(f"forceMag = {forceMag}")
    print(f"f = {f}")
    print(f"Drude受到的力：{f[0]} (向左)")
    
    print("\n\n真正的问题可能是：")
    print("-"*40)
    print("这不是库仑力计算的bug，而是整体力平衡的问题！")
    print("")
    print("在SCF迭代中：")
    print("1. Drude受到来自External的库仑力")
    print("2. Drude还受到弹簧力")
    print("3. 可能还有其他力？")
    print("")
    print("也许问题在于：")
    print("- 力的初始化")
    print("- 力的累加方式")
    print("- 或者某些力被重复计算")

def propose_debugging():
    """提出调试方案"""
    print("\n\n调试方案")
    print("="*60)
    
    print("1. 在calculateForces中打印每种力的贡献：")
    print("   - 弹簧力")
    print("   - 屏蔽库仑力")
    print("   - 普通库仑力")
    print("")
    print("2. 检查力的初始化：")
    print("   - forces是否正确清零？")
    print("   - 是否有力被重复添加？")
    print("")
    print("3. 对比FBP的力计算：")
    print("   - FBP直接计算外力")
    print("   - SCF计算所有力")
    print("   - 可能有差异")
    print("")
    print("4. 创建最小测试：")
    print("   - 只有一个Drude和一个外部电荷")
    print("   - 逐步追踪力的计算")

if __name__ == "__main__":
    analyze_physics_conventions()
    check_our_implementation()
    propose_debugging()