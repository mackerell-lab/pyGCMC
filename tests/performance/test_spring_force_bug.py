#!/usr/bin/env python3
"""
发现SCF中弹簧力的bug
"""

import numpy as np

def analyze_spring_force_bug():
    """分析弹簧力的bug"""
    print("SCF弹簧力Bug分析")
    print("="*60)
    
    print("\n查看DrudeForce.cpp第335-337行的代码：")
    print("```cpp")
    print("// Force on Drude: F = -k * (r_drude - r_parent) = -k * delta")
    print("Vec3 f = delta * drude.kIsotropic;")
    print("forces[p] -= f;  // Force on Drude (restoring force)")
    print("```")
    
    print("\n问题分析：")
    print("-"*40)
    print("1. 注释说：F = -k * delta")
    print("2. 但代码写的是：f = delta * k (没有负号！)")
    print("3. 然后：forces[p] -= f")
    print("4. 结果：forces[p] = forces[p] - k*delta = forces[p] - k*(r_d - r_p)")
    print("")
    print("这相当于：F_spring = -k*(r_d - r_p) ✓ 看起来正确")
    
    print("\n但是等等！让我们仔细分析力的更新：")
    print("-"*40)
    
    # 数值例子
    k = 138935.0
    
    print("\n情况1：Drude在Parent右边（delta > 0）")
    drude_x = 0.001
    parent_x = 0.0
    delta = drude_x - parent_x  # 0.001
    
    print(f"  delta = {delta} nm")
    print(f"  f = delta * k = {delta * k} kJ/(mol·nm)")
    print(f"  forces[drude] -= f  →  forces[drude] -= {delta * k}")
    print(f"  实际弹簧力 = -{delta * k} = {-delta * k} kJ/(mol·nm)")
    print(f"  方向：向左 ✓ 正确（拉回parent）")
    
    print("\n情况2：如果代码中有其他地方初始化了forces...")
    print("让我们看看calculateCoulombEnergy是如何处理力的")
    
    print("\n\n关键洞察：")
    print("="*60)
    print("SCF的弹簧力计算本身可能是正确的！")
    print("问题可能在于：")
    print("1. 静电力的计算")
    print("2. 力的初始化")
    print("3. 力的组合方式")
    print("4. 或者在calculateCoulombEnergy中的符号约定")

def compare_force_conventions():
    """比较不同的力约定"""
    print("\n\n力约定对比")
    print("="*60)
    
    print("物理约定：")
    print("- 弹簧力：F = -k*(r_d - r_p)")
    print("- 如果Drude在右(r_d > r_p)，力向左(F < 0)")
    print("")
    
    print("SCF代码：")
    print("```cpp")
    print("Vec3 f = delta * k;      // f = k*(r_d - r_p)")
    print("forces[p] -= f;          // forces = forces - f")
    print("```")
    print("如果forces初始为0，则forces[p] = -k*(r_d - r_p) ✓")
    print("")
    
    print("FBP代码：")
    print("```cpp")
    print("Vec3 displacement = totalForces[i] * (1.0 / k);")
    print("```")
    print("这里totalForces应该是外力，位移 = 外力/k")
    
    print("\n\n真正的问题可能在于calculateCoulombEnergy！")

def test_force_calculation_order():
    """测试力计算的顺序"""
    print("\n\n力计算顺序分析")
    print("="*60)
    
    print("SCF中calculateForces的顺序：")
    print("1. calculateHarmonicEnergy - 计算弹簧力")
    print("2. calculateScreenedCoulombEnergy - 计算屏蔽的库仑力")
    print("3. calculateCoulombEnergy - 计算所有库仑力")
    print("")
    print("问题：calculateCoulombEnergy可能包含了不该包含的相互作用！")
    print("或者符号约定不一致！")

if __name__ == "__main__":
    analyze_spring_force_bug()
    compare_force_conventions()
    test_force_calculation_order()