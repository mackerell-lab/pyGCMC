#!/usr/bin/env python3
"""
测试力的方向和平衡
"""

import numpy as np

def analyze_force_directions():
    """分析力的方向"""
    print("力方向分析")
    print("="*60)
    
    # 系统设置
    print("系统配置:")
    print("  Parent原子O: 位置(0,0,0), 电荷+1.71636")
    print("  Drude粒子D: 初始位置(0,0,0), 电荷-1.71636")
    print("  外部原子E: 位置(0.5,0,0), 电荷+1.71636")
    print("")
    
    # 参数
    q_parent = 1.71636
    q_drude = -1.71636
    q_external = 1.71636
    k_spring = 418400.0
    k_elec = 138.935
    
    # 分析两种位移情况
    displacements = [-0.003853, 0.003913]  # SCF和FBP的结果，单位nm
    names = ["SCF", "FBP"]
    
    for disp, name in zip(displacements, names):
        print(f"\n{name}结果分析 (位移={disp*1000:.3f} pm):")
        print("-"*50)
        
        # Drude最终位置
        drude_x = 0.0 + disp
        print(f"  Drude最终位置: x = {drude_x:.6f} nm")
        
        # 1. 弹簧力分析
        print(f"\n  弹簧力:")
        # F_spring = -k * (r_drude - r_parent)
        # 如果Drude在Parent右边(disp>0)，弹簧被拉伸，力指向左(负)
        # 如果Drude在Parent左边(disp<0)，弹簧被压缩，力指向右(正)
        F_spring = -k_spring * disp
        print(f"    位移向量: Δr = {disp:.6f} nm")
        print(f"    弹簧力: F_spring = -k*Δr = {F_spring:.1f} kJ/(mol·nm)")
        print(f"    方向: {'←(向左)' if F_spring < 0 else '→(向右)'}")
        
        # 2. 电场力分析
        print(f"\n  电场力:")
        
        # 外部正电荷在Drude处产生的电场
        r_vec = np.array([drude_x - 0.5, 0, 0])  # 从外部电荷指向Drude
        r_mag = np.linalg.norm(r_vec)
        r_hat = r_vec / r_mag
        
        # 电场：正电荷产生的电场沿径向向外
        E_field = k_elec * q_external / r_mag**2 * r_hat
        print(f"    从E到D的矢量: r = ({r_vec[0]:.6f}, 0, 0) nm")
        print(f"    距离: |r| = {r_mag:.6f} nm")
        print(f"    电场: E = {E_field[0]:.1f} kJ/(mol·nm·e) {'←' if E_field[0] < 0 else '→'}")
        
        # Drude受到的力
        F_electric = q_drude * E_field[0]
        print(f"    电力: F = q*E = {q_drude:.3f} * {E_field[0]:.1f}")
        print(f"         = {F_electric:.1f} kJ/(mol·nm) {'←' if F_electric < 0 else '→'}")
        
        # 3. 力平衡检查
        print(f"\n  力平衡:")
        F_total = F_spring + F_electric
        print(f"    弹簧力: {F_spring:.1f} kJ/(mol·nm)")
        print(f"    电力: {F_electric:.1f} kJ/(mol·nm)")
        print(f"    总力: {F_total:.1f} kJ/(mol·nm)")
        print(f"    {'✓ 力平衡' if abs(F_total) < 1.0 else '✗ 力不平衡'}")
        
        # 4. 物理合理性
        print(f"\n  物理合理性:")
        if disp < 0:
            print(f"    Drude向左移动(x减小)")
            print(f"    这意味着被右侧的正电荷排斥")
            print(f"    对于负电荷，这是{'✗ 不合理' if True else '✓ 合理'}的")
        else:
            print(f"    Drude向右移动(x增大)")
            print(f"    这意味着被右侧的正电荷吸引")
            print(f"    对于负电荷，这是{'✓ 合理' if True else '✗ 不合理'}的")
    
    # 理论分析
    print("\n\n理论分析:")
    print("="*60)
    print("在平衡位置，应该有 F_spring + F_electric = 0")
    print("")
    print("设Drude位移为d(向右为正):")
    print("- 弹簧力: F_spring = -k*d")
    print("- 电场力: F_electric = q_drude * E")
    print("- 其中E取决于Drude位置")
    print("")
    print("如果外部是正电荷，Drude是负电荷:")
    print("- 负电荷应该被吸引，向右移动(d>0)")
    print("- 弹簧力向左(-k*d<0)")
    print("- 电力向右(应该>0)")
    print("- 两者平衡")
    
    # 分析为什么结果不同
    print("\n\n结果差异分析:")
    print("="*60)
    print("SCF: 给出负位移，表明算法可能在最小化不同的能量函数")
    print("FBP: 给出正位移，符合物理直觉")
    print("")
    print("可能的原因:")
    print("1. SCF可能包含了Parent-Drude相互作用(不应该包含)")
    print("2. 力的计算或符号定义可能不一致")
    print("3. 边界条件或约束的处理不同")

if __name__ == "__main__":
    analyze_force_directions()